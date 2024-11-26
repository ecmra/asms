
use crate::utils;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;

use crate::matrix;

pub struct Stats {
    pub contig: String,
    pub gpos:Vec<i64>,
    pub minpos: i64,
    pub maxpos: i64,
    pub n: usize,
    pub d: usize,
    pub pop0: u16,
    pub pop1: u16,
    pub median_abs_diff_notna_0_1: f32,
    pub median_meth0_notna: f32,
    pub median_meth1_notna: f32,
    pub n00:Vec<u16>,
    pub n01:Vec<u16> ,
    pub n10:Vec<u16> ,
    pub n11:Vec<u16> ,
    pub m0:Vec<f32> ,
    pub m1:Vec<f32>
}

/// computes separate methylation per position for clusters 0,1
///
/// also computes methylation per position in total (disregarding the clusters).
/// prints m.gpos[j], n00, n10, m0, n01, n11, m1
pub fn compute_stats(m: &matrix::Matrix,
		     clusters: &Vec<u8>) -> Stats {
    // total depth and number of converted reads per cluster per position
    let mut cln: Vec<u16> = vec![0; 2 * m.d];
    let mut clx: Vec<u16> = vec![0; 2 * m.d];
    /*
    == 4 vector of counts ==
    n00: # of 0s in cluster 0
    n01: # of 0s in cluster 1
    n10: # of 1s in cluster 0
    n11: # of 1s in cluster 1
     */
    let mut pop0:u16 = 0;
    let mut pop1:u16 = 0;
    let mut n00:Vec<u16> = vec![0;m.gpos.len()];
    let mut n01:Vec<u16> = vec![0;m.gpos.len()];
    let mut n10:Vec<u16> = vec![0;m.gpos.len()];
    let mut n11:Vec<u16> = vec![0;m.gpos.len()];
    let mut m0:Vec<f32>  = vec![0.0;m.gpos.len()];
    let mut m1:Vec<f32>  = vec![0.0;m.gpos.len()];
    for cl in clusters {
	if *cl == 0 { pop0 += 1 }
	if *cl == 1 { pop1 += 1 }
    }
    for l in 0..m.idxreads.len() {
        let cl = clusters[m.idxreads[l]] as usize;
	let j = m.idxgpos[l];
        cln[cl * m.d + j] += 1;
        clx[cl * m.d + j] += m.m[m.idxreads[l] * m.d + j] as u16;
    }
    let minpos = m.gpos[0];
    let maxpos = m.gpos[m.d - 1];
    let mut absidx: usize = 0;
    let mut absdiff: Vec<f32> = vec![0.0; m.d];
    let ( mut meth0, mut meth1 ) = ( vec![0.0f32; m.d], vec![0.0f32; m.d] );
    let ( mut diff, median ): (f32, f32);
    for j in 0..m.d {
        n00[j] = cln[j] - clx[j];
        n10[j] = clx[j];
        m0[j] = n10[j] as f32 / (n00[j] + n10[j]) as f32;
        n01[j] = cln[m.d + j] - clx[m.d + j];
        n11[j] = clx[m.d + j];
        m1[j] = n11[j] as f32 / (n01[j] + n11[j]) as f32;
        diff = f32::abs(m1 [j]- m0[j]);
        if !diff.is_nan() {
            absdiff[absidx] = diff;
            meth0[absidx] = m0[j];
            meth1[absidx] = m1[j];
            absidx += 1;
        }
    }
    if 0 == absidx {
        median = f32::NAN;
    } else {
        absdiff = absdiff[..absidx].to_vec();
        meth0 = meth0[..absidx].to_vec();
        meth1 = meth1[..absidx].to_vec();
        absdiff.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if absidx % 2 == 1 {
            median = absdiff[absidx / 2];
        } else {
            median = 0.5 * (absdiff[absidx / 2 - 1] + absdiff[absidx / 2]);
        }
    }
    let medmeth0 = utils::median(&mut meth0);
    let medmeth1 = utils::median(&mut meth1);
    return Stats {
        contig: String::from(&m.contig),
	gpos:m.gpos.to_vec(),
	minpos: minpos,
        maxpos: maxpos,
        n: m.n,
        d: m.d,
        pop0: pop0,
        pop1: pop1,
        median_abs_diff_notna_0_1: median,
        median_meth0_notna: medmeth0,
        median_meth1_notna: medmeth1,
	n00:n00,
	n01:n01,
	n10:n10,
	n11:n11,
	m0:m0,
	m1:m1
    };
}

/// print separate methylation per position for clusters 0,1
///
/// print gpos[j], n00, n10, m0, n01, n11, m1
pub fn print_stats<W: std::io::Write>( s:&Stats, mut handle: BufWriter<W>)  {
    utils::write_out(&mut handle,
		     &format!("#@N:{}\n#@D:{}\n#@K:2\n", s.n, s.d));
    utils::write_out(&mut handle,
		     &format!("#@CONTIG:{}\n", s.contig));
    utils::write_out( &mut handle,
        &format!("#@MINPOS:{}\n#@MAXPOS:{}\n", s.minpos, s.maxpos),
    );
    utils::write_out( &mut handle,
        &format!("#@POP0:{}\n#@POP1:{}\n", s.pop0, s.pop1),
    );
    for j in 0..s.d {
        utils::write_out(
            &mut handle,
            &format!(
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}\n",
                s.gpos[j],
		s.n00[j], s.n10[j], s.m0[j],
		s.n01[j], s.n11[j], s.m1[j]
            ),
        );
    }
    utils::write_out(
        &mut handle,
        &format!("#@MEDIAN_ABS_DIFF_NOTNA_0_1:{:.3}\n",
		 s.median_abs_diff_notna_0_1),
    );
    utils::write_out(
        &mut handle,
        &format!("#@MEDIAN_METH0_NOTNA:{:.3}\n",
	s.median_meth0_notna),
    );
    utils::write_out(
        &mut handle,
        &format!("#@MEDIAN_METH1_NOTNA:{:.3}\n",
	s.median_meth1_notna),
    );
}

/// create summary stats file
///
/// read the comment fields #@ in the cluster files
/// and collect them all in a summary.
pub fn print_summary<W: std::io::Write>(bedfn: &str, outdir: &str, mut handle: BufWriter<W>) {
    let cols = [
        "CONTIG",
        "MINPOS",
        "MAXPOS",
        "N",
        "K",
        "D",
        "POP0",
        "POP1",
        "MEDIAN_ABS_DIFF_NOTNA_0_1",
        "MEDIAN_METH0_NOTNA",
        "MEDIAN_METH1_NOTNA",
    ];
    utils::write_out(&mut handle, &format!("#"));
    for k in 0..cols.len() - 1 {
        utils::write_out(&mut handle, &format!("{}\t", cols[k]));
    }
    utils::write_out(&mut handle, &format!("{}\t", cols[cols.len() - 1]));
    utils::write_out(&mut handle, &format!("filename\n"));
    let bufreader = BufReader::new(File::open(bedfn).unwrap());
    for line in bufreader.lines() {
        let line = line.unwrap();
        if line.starts_with("#") {
            continue;
        }
        let fields: Vec<&str> = line.trim_end().split_whitespace().collect::<Vec<_>>();
        let contig = fields[0];
        let start: i64 = fields[1].parse::<i64>().unwrap();
        let end: i64 = fields[2].parse::<i64>().unwrap();
        let statsfn = &format!("{}/stats/stats-{}-{}-{}.txt", outdir, contig, start, end);
        let check = Path::try_exists(Path::new(statsfn)).unwrap();
        if !check {
            eprintln!("skipping file {} (does not exist)", statsfn);
            continue;
        }
        let meta = map_of_stats(statsfn);
        for k in 0..cols.len() {
            utils::write_out(&mut handle, &format!("{}\t", meta[cols[k]]));
        }
        utils::write_out(&mut handle, &format!("{}\n", statsfn));
    }
}

/// import metadata from a stats file into an hashmap
///
/// keys of the hashmap are in the cols variable above
fn map_of_stats(statsfn: &str) -> FxHashMap<String, String> {
    let mut s: FxHashMap<String, String> = FxHashMap::default();
    //eprintln!("statsfn={}", statsfn);
    let bufreader = BufReader::new(File::open(statsfn).expect(&format!("Can't find {}", statsfn)));
    for line in bufreader.lines() {
        let line = line.unwrap();
        if line[0..2] != *"#@" {
            continue;
        }
        let idx = line.find(':').unwrap();
        let key = line[2..idx].to_string();
        s.insert(key, line[idx + 1..].to_string());
    }
    s
}

#[cfg(test)]
pub mod test {
    #[test]
    fn compute_stats_works() {
	// m is 4 X 5
	let m = crate::matrix::test::make_test_matrix();
	let cl:Vec<u8>=vec![1, 0, 1, 0];
	let s = super::compute_stats(&m, &cl);
	eprintln!("{}\t{}\t{}\t{}", s.n, s.d, s.pop0, s.pop1);
	assert!((s.pop0 == 2) && (s.pop1 == 2))
    }
}

/*
info collected by the summary of all the statistics:
cols=["MINPOS", "MAXPOS", "N", "K" ,"D", "POP0" ,"POP1" ,"MEDIAN_ABS_DIFF_NOTNA_0_1", "CONTIG"]
*/

/* output of cvlr-cluster (#@IT and LL excluded)
#@N:10
#@D:12
#@K:2
#@DF:25
#@INITPI:0.400  0.600
#@GPOS:0        1       2       3       4       5       6       7       8       9       10      11
#@INITMU0:0.100 0.100   0.100   0.100   0.100   0.100   0.900   0.100   0.900   0.100   0.100   0.100
#@INITMU1:0.100 0.100   0.100   0.100   0.100   0.100   0.100   0.100   0.100   0.100   0.100   0.100
#@LL:-60.6      BIC:-88.2       AIC:-85.6
#@PI:0.500      0.500
#@MU0:0.600     0.400   0.600   0.400   1.000   0.400   0.600   0.000   1.000   0.800   0.600   0.800
#@MU1:0.800     0.200   0.000   0.200   0.000   0.200   0.600   0.600   0.200   0.400   0.200   0.000
#@POP0:10
#@POP1:5
r0      1       0       1
*/

/// stores meth info, per read and per genomic position
///
/// n: number of reads (rows)
/// d: number of CpGs (columns)
/// m: methylation state (-1:unknown,/0/1)
/// gpos: Vec of d genomic coordinates, one per column
/// rnames: Vec of n read names, one per row
pub struct Matrix {
    pub contig: String,
    pub n: usize,
    pub d: usize,
    pub m: Vec<f32>,
    pub gpos: Vec<i64>,
    pub rnames: Vec<String>,
    /* 
    m is a sparse matrix (where absence is coded as -1).
    idxgpos, idxreads store (in parallel) the columns, rows
    of the stored values (different from -1). 
    indexing is row-wise similar to C,
    so position (i,j) maps to index i*m.d+j
     */
    pub idxgpos: Vec<usize>,
    pub idxreads: Vec<usize>,
}
use crate::commands::cluster;
use fxhash::FxHashMap;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Record;
use rust_htslib::htslib::BAM_FDUP;
use rust_htslib::htslib::BAM_FQCFAIL;
use rust_htslib::htslib::BAM_FSECONDARY;
use rust_htslib::htslib::BAM_FSUPPLEMENTARY;
use rust_htslib::htslib::BAM_FUNMAP;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufWriter;

/// create a matrix containing the methylation status
///
/// of each CpG on each considered read.
/// rows are reads, columns are genomic positions.
/// the matrix is written in a file inside the `outdir` directory.
pub fn meth_of_bam(
    bamreader: &mut IndexedReader,
    r: cluster::Region,
    outdir: &str,
    refseq: &[u8],
) -> Matrix {
    let mut discarded: u32 = 0;
    let contig = r.contig;
    let start = r.start;
    let end = r.end;
    // from the docs:
    // "start and stop are zero-based. start is inclusive, stop is exclusive."
    // https://docs.rs/rust-htslib/0.47.0/rust_htslib/bam/struct.IndexedReader.html#method.fetch
    bamreader.fetch((contig, start, end)).unwrap();
    let mut mod_in_read: u32;
    let mut mod_count: u32 = 0;
    let mut rcount: usize = 0;
    let mut rnames = FxHashMap::default();
    let mut cpg_gpos: Vec<i64> = Vec::with_capacity(1000);
    let mut cpg_read: Vec<usize> = Vec::with_capacity(1000);
    let mut cpg_state: Vec<i8> = Vec::with_capacity(1000);
    let mut cpg_gpos_set = BTreeSet::new();
    let mut recordbuf = Record::new();
    while let Some(result) = bamreader.read(&mut recordbuf) {
        let record = match result {
            Ok(_) => &recordbuf,
            Err(_) => panic!("BAM parsing failed"),
        };
        let isrev = record.is_reverse();
        let rname = std::str::from_utf8(record.qname()).unwrap();
        if skip(&record) {
            discarded += 1;
            continue;
        }
        let ap = record.aligned_pairs();
        // fill a map rpos => gpos.
        // gpos is 0-based
        let mut rpos2gpos = FxHashMap::default();
        for p2 in ap {
            rpos2gpos.insert(p2[0], p2[1]);
        }
        let mods = match record.basemods_iter() {
            Ok(mods) => mods,
            Err(_) => panic!("Failed to parse MM/ML aux tags in {}\n", rname),
        };
        // iterate over the modifications in this record
        mod_in_read = 0;
        for res in mods {
            if let Ok((position, m)) = res {
                let pos = position as i64;
                if (rpos2gpos.contains_key(&pos)) && ('m' == m.modified_base as u8 as char) {
                    let gpos = if isrev {
                        rpos2gpos[&pos]
                    } else {
                        rpos2gpos[&pos] + 1
                    };
                    //eprintln!("rpos:{} gpos:{}", pos, gpos);
                    if !validate_cgmatch(gpos, start, end, refseq) {
                        continue;
                    }
                    cpg_gpos.push(gpos);
                    cpg_read.push(rcount);
                    cpg_state.push(binarize(m.qual.try_into().unwrap(), 127, 127));
                    cpg_gpos_set.insert(gpos);
                    mod_in_read += 1;
                    mod_count += 1
                }
            }
        }
        if mod_in_read > 0 {
            rnames.insert(rcount, String::from(rname));
            rcount += 1;
        }
    }
    let outfn = &format!(
        "{}/matrices/matrix-{}-{}-{}.txt",
        outdir, contig, start, end
    );
    let handle = BufWriter::new(File::create(outfn).unwrap());
    let m = print_matrix(
        &cpg_gpos,
        &cpg_gpos_set,
        &cpg_read,
        &cpg_state,
        &rnames,
        discarded,
        mod_count,
        r,
        handle,
    );
    return m;
}

/// skip dubious reads
///
/// based on flags, presence of MM, quality must be >= 10
fn skip(record: &Record) -> bool {
    let flags = record.flags();
    let rname = std::str::from_utf8(record.qname()).unwrap();
    if flags & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP) as u16
        != 0
    {
        eprintln!("skipping {} [flag:{}]", rname, flags);
        return true;
    }
    if !(record.aux(b"MM").is_ok() || record.aux(b"Mm").is_ok()) {
        eprintln!("skipping read {} (no MM/ML aux tag)", rname);
        return true;
    }
    if record.mapq() < 10 {
        eprintln!("skipping {} [qual:{}]", rname, record.mapq());
        return true;
    }
    return false;
}

/// associate -1,0,1 to a meth value between 0 and 255
///
/// meth in (lb,ub] is called as ambiguous/unknown (-1).
fn binarize(meth: u8, lb: u8, ub: u8) -> i8 {
    if meth <= lb {
        return 0;
    }
    if (meth > lb) & (meth <= ub) {
        return -1;
    }
    if meth > ub {
        return 1;
    }
    panic!("error: illegal meth score:{}\n", meth);
}

/// check that a given gpos corresponds to a CG on the reference
fn validate_cgmatch(gpos: i64, start: i64, end: i64, refseq: &[u8]) -> bool {
    let check1: bool = gpos >= start && gpos <= end;
    let check2: bool;
    // gpos is 0-based
    if (gpos > 0) && (gpos < refseq.len() as i64) {
        check2 = b'C' == refseq[(gpos - 1) as usize] && b'G' == refseq[gpos as usize];
    } else {
        check2 = false;
    }
    return check1 && check2;
}
use crate::utils;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Read;

/// print matrix and return a matrix struct
///
/// only combinations of (read, gpos) where the methylation
/// status is known are printed. The rest of the matrix is implicitly
/// assumed to be -1.
fn print_matrix<W: std::io::Write>(
    cpg_gpos: &Vec<i64>,
    cpg_gpos_set: &BTreeSet<i64>,
    cpg_read: &Vec<usize>,
    cpg_state: &Vec<i8>,
    rnames: &FxHashMap<usize, String>,
    discarded: u32,
    mod_count: u32,
    region: cluster::Region,
    mut handle: BufWriter<W>,
) -> Matrix {
    let vecgpos: Vec<i64> = cpg_gpos_set.clone().into_iter().collect();
    let mut vecrnames: Vec<String> = vec!["".to_string(); rnames.len()];
    let mut idxgpos: Vec<usize> = vec![0; cpg_gpos.len()];
    for i in 0..rnames.len() {
        vecrnames[i] = rnames[&i].clone();
    }
    let contig = region.contig;
    utils::write_out(&mut handle, &format!("#@N:{}\n", rnames.len()));
    utils::write_out(&mut handle, &format!("#@D:{}\n", cpg_gpos_set.len()));
    utils::write_out(&mut handle, &format!("#@CONTIG:{}\n", contig));
    utils::write_out(
        &mut handle,
        &format!("#@REGION:{}:{}-{}\n", contig, region.start, region.end),
    );
    utils::write_out(&mut handle, &format!("#@CGCOUNT:{}\n", mod_count));
    utils::write_out(&mut handle, &format!("#@DISCARDED_READS:{}\n", discarded));
    if vecgpos.len() == 0 {
        return Matrix {
            contig: contig.to_string(),
            n: 0,
            d: 0,
            m: vec![0.0],
            gpos: vec![],
            rnames: vec![],
            idxgpos: vec![],
            idxreads: vec![],
        };
    }
    let n = vecrnames.len();
    let d = vecgpos.len();
    /*
    m is sparse. -1 indicates measured combination
    (read, gpos)
    */
    let mut m: Vec<f32> = vec![-1.0; n * d];
    let minpos = vecgpos[0];
    let maxpos = vecgpos[vecgpos.len() - 1];
    utils::write_out(&mut handle, &format!("#@MINPOS:{}\n", minpos));
    utils::write_out(&mut handle, &format!("#@MAXPOS:{}\n", maxpos));
    utils::write_out(&mut handle, &format!("#@SPAN:{}\n", maxpos - minpos + 1));
    for i in 0..cpg_read.len() {
        idxgpos[i] = vecgpos.iter().position(|&x| x == cpg_gpos[i]).unwrap() as usize;
        utils::write_out(
            &mut handle,
            &format!(
                "{}\t{}\t{}\t{}\t{}\n",
                cpg_read[i], rnames[&cpg_read[i]], idxgpos[i], cpg_gpos[i], cpg_state[i]
            ),
        );
        m[cpg_read[i] * d + idxgpos[i]] = cpg_state[i] as f32;
    }
    return Matrix {
        contig: contig.to_string(),
        n: n,
        d: d,
        m: m,
        gpos: vecgpos.clone(),
        rnames: vecrnames.clone(),
        idxgpos: idxgpos,
        idxreads: cpg_read.clone(),
    };
}

#[allow(dead_code)]
fn make_matrix(contig: &String, gpos: &Vec<i64>, rnames: &Vec<String>, m: &Vec<f32>) -> Matrix {
    let n = rnames.len();
    let d = gpos.len();
    assert!(n * d == m.len());
    let matrix: Matrix = Matrix {
        contig: contig.to_string(),
        n: n,
        d: d,
        m: m.to_vec(),
        gpos: gpos.to_vec(),
        rnames: rnames.to_vec(),
        idxgpos: vec![],
        idxreads: vec![],
    };
    return matrix;
}

use std::io::BufRead;
use std::io::BufReader;

 
pub fn matrix_of_file(fname: &str) -> Matrix {
    let bufreader = BufReader::new(File::open(fname).unwrap());
    let mut fields: Vec<&str>;
    let mut idxreads:Vec<usize> = vec![];
    let mut idxgpos:Vec<usize> = vec![];
    let mut rnames:Vec<String> = vec![];
    let mut gpos:Vec<i64> = vec![];
    let md = utils::read_metadata(fname);
    let n = md["N"].parse::<usize>().unwrap();
    let d = md["D"].parse::<usize>().unwrap();
    let contig = md["CONTIG"].to_string();
    let mut m:Vec<f32> = vec![ -1.0; n*d ];
    for line in bufreader.lines(){
	let s = line.unwrap();
	if &s[0..2] == "#@" {continue;}
	fields = s.split('\t').collect();
	let idxr = fields[0].parse::<usize>().unwrap();
	idxreads.push(idxr);
	if !( rnames.contains(&fields[1].to_string())) {
	    rnames.push(fields[1].to_string());
	}
	let idxg = fields[2].parse::<usize>().unwrap();
	idxgpos.push(idxg);
	if !( gpos.contains(&fields[3].parse::<i64>().unwrap())) {
	    gpos.push(fields[3].parse::<i64>().unwrap());
	}
	m[idxr*d+idxg] = fields[4].parse::<f32>().unwrap();
    }
    let m = Matrix {
	contig:contig,
        n: n,
        d: d,
        m: m,
        gpos: gpos,
        rnames: rnames,
        idxgpos: idxgpos,
        idxreads: idxreads,
    };
    return m;
}

#[cfg(test)]
pub mod test {
    pub fn make_test_matrix() -> super::Matrix {
        let tm1: super::Matrix = super::Matrix {
            contig: "test".to_string(),
            n: 4,
            d: 5,
            m: vec![
                0.0, 1.0, 0.0, 1.0, 0.0,
		1.0, 0.0, 1.0, 0.0, 1.0,
		0.0, 1.0, 0.0, 1.0, 0.0,
		1.0, 0.0, 1.0, 0.0, 1.0,
            ],
            gpos: vec![1, 2, 3, 4, 5],
            rnames: vec![
                "r1".to_string(),
                "r2".to_string(),
                "r3".to_string(),
                "r4".to_string(),
            ],
            idxgpos: vec![0,1,2,3,4],
            idxreads: vec![0,1,2,3],
        };
        return tm1;
    }
    
    #[test]
    fn matrix_of_file_works() {
	let fname = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/matrices/matrix-chr15-100863815-100863868.txt";
	let m = super::matrix_of_file(fname);
	assert!( m.idxgpos[1] == 2 );
	assert!( m.idxreads[1] == 0 );
	assert!( m.m[0*m.d+1] == 0.0 );
    }
}

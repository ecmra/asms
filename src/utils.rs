#![allow(dead_code)]

use crate::matrix::Matrix;
use rust_htslib::tbx::{self, Read};
use std::io::BufWriter;
use std::io::Write;
use std::str::FromStr;

pub struct MethSlice {
    pub contig: String,
    pub pos: Vec<u64>,
    pub modif: Vec<u32>,
    pub canon: Vec<u32>,
    pub depth: Vec<u32>,
    pub meth: Vec<f32>,
}

use std::io::BufReader;
use std::fs::File;
use fxhash::FxHashMap;
use std::io::BufRead;
pub fn read_metadata(fname:&str) -> fxhash::FxHashMap<String, String> {
    let mut bufreader = BufReader::new(File::open(fname).unwrap());
    let mut metadata:FxHashMap<String, String> = FxHashMap::default();
    let mut fields: Vec<&str>;
    let mut line = String::new();
    loop {
	let _ = bufreader.read_line(&mut line); 
        if &line[0..2] == "#@" {
            fields = line.split(':').collect();
            fields[0] = &fields[0][2..];
            metadata.insert(fields[0].to_string(), fields[1].trim().to_string());
	    line.clear();
        } else {
	    break;
	}
    }
    return metadata;
}



pub fn write_out<W: std::io::Write>(handle: &mut BufWriter<W>, s: &str) {
    handle.write_all(s.as_bytes()).unwrap();
}

/// load up a bedmethyl region in memory.
///
/// only the CpGs, only methylation data.
/// According to specs [https://nanoporetech.github.io/modkit/intro_bedmethyl.html]
/// column 4 contains a
/// "single letter code for modified base and motif when more than one motif is used"
/// here we look for 'm' in fields[3] (4th column).
///
/// fields[5] (the 6th column) corresponds to strand.
/// Here I assume I am working with combined strands, hence
/// I want strand to be '.'. On the other hand, if the strands are not combined
/// this is not reason enough to stop
/// processing a bedmethyl file. If the strand is not '.' then only positions on the
/// '+' strand are considered, with a warning.
///
/// the last version of modkit (0.3.0, f9320a) changes the output format.
/// previous format has tab separated columns for the first 10 fields and space thereafter.
/// example bedmethyl line:
/// chr1\t10468\t10469\tm\t15\t.\t10468\t10469\t255,0,0\t15 80.00 12 3 0 0 6 2 0
/// the new format has tabs everywhere.
/// they are both ok by the standard,
/// https://samtools.github.io/hts-specs/BEDv1.pdf
/// although only the new format follows the recommendation of using a single
/// delimiter type.
pub fn get_bedmethyl_slice(
    tbxreader: &mut tbx::Reader,
    contig: &str,
    start: u64,
    end: u64,
) -> MethSlice {
    let tid = tbxreader.tid(&contig).unwrap();
    tbxreader.fetch(tid, start, end).unwrap();
    let mut records = tbxreader.records();
    let mut pos: Vec<u64> = vec![0; 5_000_000];
    let mut modif: Vec<u32> = vec![0; 5_000_000];
    let mut canon: Vec<u32> = vec![0; 5_000_000];
    let mut r: Vec<u8>;
    let mut s: String;
    let mut smodif: &str;
    let mut scanon: &str;
    let mut ncpg: usize = 0;
    let mut fields: Vec<&str>;
    let mut strandwarn: bool = false;
    loop {
        r = match records.next() {
            Some(r) =>
            // here it's a Result<T,E> -- unpack --
            {
                match r {
                    Ok(r) => r,
                    Err(_e) => break,
                }
            }
            None => break,
        };
        s = String::from_utf8(r).unwrap();
        fields = s.split('\t').collect();
        if fields[3] != "m" {
            continue;
        }
        if !(fields[5] == "." || fields[5] == "+") {
            strandwarn = true;
            continue;
        }
        pos[ncpg] = u64::from_str(&fields[1]).unwrap();
        if fields[9].contains(' ') {
            smodif = fields[9].split_whitespace().collect::<Vec<&str>>()[2];
            scanon = fields[9].split_whitespace().collect::<Vec<&str>>()[3];
        } else {
            smodif = fields[11];
            scanon = fields[12];
        }
        modif[ncpg] = u32::from_str(smodif).unwrap();
        canon[ncpg] = u32::from_str(scanon).unwrap();
        ncpg += 1;
    }
    if strandwarn {
        eprintln!("Warning: positions with negative strand have been skipped");
    }
    let modif = modif[0..ncpg].to_vec();
    let canon = canon[0..ncpg].to_vec();
    let depth: Vec<u32> = modif.iter().zip(canon.iter()).map(|(m, c)| m + c).collect();
    let meth: Vec<f32> = modif
        .iter()
        .zip(depth.iter())
        .map(|(m, d)| *m as f32 / *d as f32)
        .collect();
    MethSlice {
        contig: String::from(contig),
        pos: pos[0..ncpg].to_vec(),
        modif: modif,
        canon: canon,
        depth: depth,
        meth: meth,
    }
}

/// parse chr[:start-end] into a triplet
pub fn parse_region(region: &str) -> (&str, i64, i64) {
    let fields: Vec<&str> = region.split(":").collect();
    let contig = &fields[0];
    if 1 == fields.len() {
        return (contig, 0, i64::MAX);
    }
    let fields: Vec<&str> = fields[1].split("-").collect();
    let start = fields[0].parse::<i64>().unwrap();
    let end = fields[1].parse::<i64>().unwrap();
    return (contig, start, end);
}

/// median of Vec<f32>
pub fn median(v: &mut Vec<f32>) -> f32 {
    let vl = v.len();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if vl % 2 == 1 {
        return v[vl / 2];
    } else {
        return 0.5 * (v[vl / 2 - 1] + v[vl / 2]);
    }
}

/// hamming distance
///
/// input vectors are ternary: ecah position can be -1,0,1
/// positions which are -1 in either of the input vectors
/// are skipped in the computation of the distance
/// return hamming distance and number of positions considered.
pub fn hammingd(v1: &[f32], v2: &[f32]) -> (usize, usize) {
    let d = v1.len();
    let mut hd: usize = 0;
    let mut n: usize = 0;
    for i in 0..d {
        if v1[i] == -1.0 || v2[i] == -1.0 {
            continue;
        }
        if v1[i] != v2[i] {
            hd += 1;
            n += 1;
            continue;
        }
        n += 1;
    }
    return (hd, n);
}

/// compute hamming distance between rows of matrix
///
/// r1, r2 row indexes
/// return hamming distance and number of positions considered
pub fn hammingd_rows(m: &Matrix, r1: usize, r2: usize) -> (usize, usize) {
    let v1 = &m.m[r1 * m.d..m.d * (r1 + 1)];
    let v2 = &m.m[r2 * m.d..m.d * (r2 + 1)];
    return hammingd(v1, v2);
}

/// compute average intra and inter distance between clusters of rows
///
/// return intra, inter
pub fn intra_inter_hammingd(m: &Matrix, cl: &Vec<u8>) -> (f32, f32) {
    let mut nintra: usize = 0;
    let mut dintra: f32 = 0.0;
    let mut ninter: usize = 0;
    let mut dinter: f32 = 0.0;
    for i in 0..m.n {
        for j in i..m.n {
            let d = hammingd_rows(&m, i, j).0;
            if cl[i] == cl[j] {
                nintra += 1;
                dintra += d as f32;
            } else {
                ninter += 1;
                dinter += d as f32;
            }
        }
    }
    return (dintra / nintra as f32, dinter / ninter as f32);
}

/// compute mean silhouette score
///
/// at the moment stats are computed inside cluster.rs
/// hence this fuction should be called from there.
pub fn mean_silhouette_score(m: &Matrix, cl: &Vec<u8>) -> f32 {
    let (avgintra, avginter) = intra_inter_hammingd(m, cl);
    return (avginter - avgintra) / f32::max(avginter, avgintra);
}

#[test]
fn hammingd_test() {
    let v1 = vec![1.0, 0.0, 1.0];
    let v2 = vec![0.0, 1.0, 1.0];
    let res = hammingd(&v1, &v2);
    assert!(res.0 == 2 && res.1 == 3);
    let v1 = vec![-1.0, -1.0, 0.0];
    let v2 = vec![0.0, -1.0, 1.0];
    let res = hammingd(&v1, &v2);
    assert!(res.0 == 1 && res.1 == 1);
    // explicit slices
    let a = vec![0.0, 1.0, 0.0, 1.0, 0.0, 0.0];
    let a1 = &a[0..3];
    let a2 = &a[3..];
    let res = hammingd(a1, a2);
    assert!(res.0 == 2 && res.1 == 3);
}


#[cfg(test)]
pub mod test {
#[test]
    fn read_metadata_works() {
	let fname = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/matrices/matrix-chr15-100863815-100863868.txt";
	let md = super::read_metadata(fname);
	assert!(md["N"] == "96".to_string());
	assert!(md["D"] == "3".to_string());
    }
}

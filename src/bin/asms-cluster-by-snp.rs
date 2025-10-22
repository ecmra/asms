
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;

use rust_htslib::bam::Record;

use fxhash::FxHasher;
use std::hash::BuildHasherDefault;
use std::collections::HashMap;
use rust_htslib::bam::ext::BamRecordExtensions;


use rust_htslib::htslib::BAM_FDUP;
use rust_htslib::htslib::BAM_FQCFAIL;
use rust_htslib::htslib::BAM_FSECONDARY;
use rust_htslib::htslib::BAM_FSUPPLEMENTARY;
use rust_htslib::htslib::BAM_FUNMAP;

use asms::parsemod;

//use numutils::ch14;

/*

test with

target/release/asms-cluster-by-snp /home/ecmra/data/ashkenazi-trio-chr15/son_HG002.chr15.bam  ../asms-tests/out/cluster-by-snp-test/variants.txt

*/

fn main() {
    /*
    
    SNPs (contig:gcoord:ref:alt)
     */
    let args:Vec<String> = env::args().collect();
    let bamfn= &args[1];
    let varfn = &args[2];
    let stdin= std::io::stdin();
    let varf:Box<dyn BufRead>;
    if "-" == varfn {
    	varf = Box::new(stdin.lock());
    } else {
    	varf = Box::new(BufReader::new(File::open(varfn).unwrap()));
    }
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    let mut recordbuf = Record::new();
    let mut meth:Vec<f32>;
    let mut lc=0;
    for line in varf.lines() {
    	let line = line.unwrap();
	lc+=1;
    	eprintln!("{lc}:{line}");
    	let fields:Vec<&str> = line.split(':').collect();
	let contig = fields[0];
	let gpos1:i64 = fields[1].parse::<i64>().unwrap();
	let start:i64 = gpos1 - 500 ;
	let end:i64 = gpos1 + 500 ;
	let refnuc = fields[2].bytes().nth(0).unwrap();
	let altnuc = fields[3].bytes().nth(0).unwrap();
	bamreader.fetch((contig, start, end)).unwrap();
	while let Some(result) = bamreader.read(&mut recordbuf) {
            let record = match result {
		Ok(_) => &recordbuf,
		Err(_) => panic!("BAM parsing failed"),
            };
	    let flags = record.flags();
	    let qflags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP; 
	    if flags & qflags as u16 != 0 { continue }
	    if record.mapq() < 10 { continue }
	    let rname = std::str::from_utf8(record.qname()).unwrap();
	    let seq = record.seq();
	    let ap =  record.aligned_pairs();
	    let hasher = BuildHasherDefault::<FxHasher>::default();
	    let mut rpos0_of_gpos0:HashMap<i64,i64,BuildHasherDefault<FxHasher> > = HashMap::with_capacity_and_hasher(10000, hasher);
	    let mut minrpos0:i64 = i64::MAX;
	    let mut maxrpos0:i64 = 0;
	    // very expensive, search ap instead
	    for [rpos0, gpos0] in ap {
		let tmp1 = gpos0 + 1;
		if tmp1 >= start && tmp1 <= end {
		    if rpos0 < minrpos0 {minrpos0 = rpos0}
		    if rpos0 > maxrpos0 {maxrpos0 = rpos0}
		    rpos0_of_gpos0.insert(gpos0, rpos0);
		}
	    }
	    let gpos0 = gpos1-1;
	    if !rpos0_of_gpos0.contains_key(&gpos0) {continue}
	    let rpos0 = rpos0_of_gpos0[&gpos0];
	    let nuc = unsafe{ seq.decoded_base_unchecked(rpos0 as usize)};
	    let hp = match nuc {
		_ if nuc == refnuc => 0u8,
		_ if nuc == altnuc => 1u8,
		_ => continue,
	    };
	    meth = vec![];
	    let mut calls = vec![];
	    parsemod::parse_mods(&record, &mut calls, seq, minrpos0, maxrpos0);
            calls.sort_by(|c1, c2| c1.0.cmp(&c2.0));
	    // calls contain meth info for all the Cs not specifically the CGs.
	    for res in calls {
		let modrpos0 = res.0;
		let m = res.1;
		// C and G could be nearby on the read, and far away on the genome
		// because of a deletion; I assume this does not happen frequently,
		// hence I don't check.
		let (nuc0, nuc1) = match record.is_reverse() {
		    false => (unsafe{ seq.decoded_base_unchecked(modrpos0 as usize) }, 
			      unsafe{ seq.decoded_base_unchecked((modrpos0 + 1) as usize) }),
		    true => (unsafe{ seq.decoded_base_unchecked((modrpos0 - 1) as usize) }, 
			      unsafe{ seq.decoded_base_unchecked(modrpos0 as usize) }),
		};
		let checkcg = (nuc0 == b'C') && (nuc1 == b'G');
		if !(b'm' == m.modified_base as u8) || !checkcg {continue}
		let mqual = m.qual as u8;
		if mqual > 127 {meth.push(1.0)}
		else {meth.push(0.0)}
	    }
	    let ncpgs = meth.len();
	    let mmeth = meth.iter().sum::<f32>()/ncpgs as f32;
	    let sign = match record.is_reverse() {
		false => '+',
		true  => '-',
	    };
	    match hp {
		0 => println!("{rname}\t0\t{mmeth}\t{ncpgs}\t{line}\t{sign}"),
		1 => println!("{rname}\t1\t{mmeth}\t{ncpgs}\t{line}\t{sign}"),
		_ => panic!("EEEEEEEEEEEEE????"),
	    }
	}
    }
}

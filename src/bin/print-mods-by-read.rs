use std::env;

use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

use rust_htslib::htslib;

use asms::parsemod;

//use std::any::type_name_of_val;

fn print_calls(rname: &str,
		    calls: &mut Vec<(i32, htslib::hts_base_mod)>) {
    
    for res in calls {
        let (rpos0, m) = res;
	let strand = match m.strand {
	    0 => '+',
	    1 => '-',
	    _ => 'X'
	};
	println!(
            "{rname}:{rpos0}:{0}:{1}:{2}:{3}",
            m.canonical_base as u8 as char,
	    m.modified_base as u8 as char,
	    m.qual, strand,
        );
    }
}

fn main() {
    let cmdargs: Vec<String> = env::args().collect();
    let bamfn = &cmdargs[1];
    let contig = &cmdargs[2];
    let start: i32 = cmdargs[3].parse().unwrap();
    let end: i32 = cmdargs[4].parse().unwrap();
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    bamreader.fetch((contig, start, end)).unwrap();
    let mut recordbuf = Record::new();
    while let Some(result) = bamreader.read(&mut recordbuf) {
        let record = match result {
            Ok(_) => &recordbuf,
            Err(_) => panic!("BAM parsing failed"),
        };
        let rname = std::str::from_utf8(record.qname()).unwrap();
        let mut calls = vec![];
	let seq = record.seq();
        parsemod::parse_mods(&record, &mut calls, seq, i32::MIN, i32::MAX);
	calls.sort_by(|c1, c2| c1.0.cmp(&c2.0));
	print_calls(rname, &mut calls);
    }
}

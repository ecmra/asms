
use std::env;

use rust_htslib::bam::record::Aux;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read;


use asms::matrix;

fn get_mm<'a>(mm:&'a str, record:&'a Record) -> &'a str {
    let mm = record.aux(mm.as_bytes());
    let mm = match mm {
	Ok(Aux::String(mm)) => mm,
	_ => ""
    };
    return mm;
}

fn main() {
    let cmdargs:Vec<String> = env::args().collect();
    let bamfn = &cmdargs[1];
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    let contig = &cmdargs[2];
    let start:i32 = cmdargs[3].parse().unwrap();
    let end:i32   = cmdargs[4].parse().unwrap();
    bamreader.fetch((contig, start, end)).unwrap();
    let mut recordbuf = Record::new();
    while let Some(result) = bamreader.read(&mut recordbuf) {
        let record = match result {
            Ok(_) => &recordbuf,
            Err(_) => panic!("BAM parsing failed"),
        };
	let skip = matrix::skip(&record);
	let isrev = record.is_reverse();
        let rname = std::str::from_utf8(record.qname()).unwrap();
	let mm = get_mm("MM", &record);
	let mm = match mm {
	    "" => get_mm("Mm", &record),
	    _ => mm
	};
	let mods = match record.basemods_iter() {
	    Ok(mods) => mods,
	    Err(_) => {
		println!("Failed to parse MM/ML aux tags in {}\n", rname);
		continue
	    }
	};
	for res in mods {
            if let Ok((position, m)) = res {
		println!("{rname}:{position}:{0}:{1}:{2}:{3}:{4}",
			 &mm[0..5],
			 m.modified_base as u8 as char,
			 m.qual, skip, isrev);
	    }
	}
    }
}

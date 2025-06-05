
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
	let rskip = match matrix::skip(&record) {
	    true => "rskip",
	    false => ""
	};
	let isrev = match record.is_reverse(){
	    true => "rev",
	    false => "norev"
	};
        let rname = std::str::from_utf8(record.qname()).unwrap();
	let mm = get_mm("MM", &record);
	let mm = match mm {
	    "" => get_mm("Mm", &record),
	    _ => mm
	};
	/* mods has all the methods of an iterator,
	   with item = Result<(i32 (pos), hts_base_mod)>
	   plus extra methods to exctract info
	   about how many and what type of
	   modifications are recorded fot the record
	   I am looking at 
	*/
	let mods = match record.basemods_iter() {
	    Ok(mods) => mods,
	    Err(_) => {
		println!("Failed to parse MM/ML aux tags in {}\n", rname);
		continue
	    }
	};
	let mcodes:&[i32] = mods.recorded();
	let mut implicit:bool = true;
	// record value of implicit for m
	for i in 0..mcodes.len(){
	    let modtype = mods.query_type(mcodes[i]).unwrap();
	    let canbase = modtype.canonical as char;
	    if (modtype.canonical as char == 'C' &&
		mcodes[i] as u8 as char == 'm'
		&& modtype.implicit == 0) { implicit = false};
	}
	// here depending on implicit I can iterate in different ways
	if implicit {
	    // this iterates over methylated positions
	    for res in mods {
		if let Ok((position, m)) = res {
		    if m.modified_base as u8 as char != 'm' { continue}
		    println!("{rname}:{position}:{0}:{1}:{2}:{3}:{4}:{5}:{6}",
			     &mm[0..5],
			     m.modified_base as u8 as char,
			     m.qual,  rskip, isrev,
			     mcodes.len(), implicit);
		}
	    }
	} else {
	    // deal with implicit calls
	    // first look at explicit calls
	    let mut calls:Vec<(i32, hts_base_mod)> = vec![];
	    for res in mods {
		if let Ok((position, m)) = res {
		    if m.modified_base as u8 as char == 'm' { calls.push(res) }
		}
	    }
	    // add implicit calls
	}
    }
}

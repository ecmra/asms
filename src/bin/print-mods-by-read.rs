use std::env;

use rust_htslib::bam::record::Aux;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

use rust_htslib::htslib;

use asms::matrix;

//use std::any::type_name_of_val;

fn get_mm<'a>(mm: &'a str, record: &'a Record) -> &'a str {
    let mm = record.aux(mm.as_bytes());
    let mm = match mm {
        Ok(Aux::String(mm)) => mm,
        _ => "",
    };
    return mm;
}

/// base modification iter
/// can give positions which correspond to a G on the read
/// when the mod is detected on the other strand.
/// for example when using C+m on a read which is
/// reverse complemented.
fn call_explicit(
    mods: rust_htslib::bam::record::BaseModificationsIter,
    calls: &mut Vec<(i32, htslib::hts_base_mod)>,
    isrev: bool,
) {
    for res in mods {
        if let Ok((position, m)) = res {
            if m.modified_base as u8 as char != 'm' {
                continue;
            }
            match (m.canonical_base as u8, isrev) {
                (b'C', false) => { calls.push((position, m)); }
                (b'C', true)  => { calls.push((position-1, m));}
                (b'G', false) => { calls.push((position-1 , m))},
		(b'G', true)  => { calls.push((position , m))}
                _ => {
                    panic!("call_explicit:panic!");
                }
            }
        }
    }
}

/// add implicit calls: scan the sequence in the 5' to 3' direction.
/// everytime there is a C, check if its position
/// (corrected for strand) is already in the mod vector. If not,
/// then push it in. Sort the resulting vector by position.
fn fill_implicit(
    record: &rust_htslib::bam::Record,
    calls: &mut Vec<(i32, htslib::hts_base_mod)>,
    canbase: char,
) {
    let isrev = record.flags() & 0x10 == 0x10;
    let explpos: Vec<i32> = calls.iter().map(|e| e.0).collect();
    let seq = record.seq().as_bytes();
    for pos in 0..record.seq_len() {
	// all the modifications correspond to Cs on the seq
	// see the function call_explicit 
	if seq[pos] == b'C' as u8 {
	    if !(explpos.contains(&(pos as i32))) {
                let mb = htslib::hts_base_mod {
                    modified_base: b'm' as i32,
                    canonical_base: seq[pos] as i32,
                    strand: match seq[pos] {
                        b'C' => 0,
                        b'G' => 1,
                        _ => panic!("fill_implicit::panic!"),
                    },
                    qual: 0,
                };
                calls.push((pos as i32, mb));
            
            }
	}
    }
}

fn sort_print_calls(
    mm: &str,
    rskip: &str,
    isrev: bool,
    mcodes: &[i32],
    implicit: bool,
    rname: &str,
    calls: &mut Vec<(i32, htslib::hts_base_mod)>,
) {
    calls.sort_by(|c1, c2| c1.0.cmp(&c2.0));
    for res in calls {
        let (position, m) = res;
        println!(
            "{rname}:{position}:{0}:{1}:{2}:{3}",
            m.canonical_base as u8 as char,
            m.modified_base as u8 as char,
            m.qual,
            m.strand,
        );
    }
}

fn main() {
    let cmdargs: Vec<String> = env::args().collect();
    let bamfn = &cmdargs[1];
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    let contig = &cmdargs[2];
    let start: i32 = cmdargs[3].parse().unwrap();
    let end: i32 = cmdargs[4].parse().unwrap();
    bamreader.fetch((contig, start, end)).unwrap();
    let mut recordbuf = Record::new();
    while let Some(result) = bamreader.read(&mut recordbuf) {
        let record = match result {
            Ok(_) => &recordbuf,
            Err(_) => panic!("BAM parsing failed"),
        };
        let rskip = match matrix::skip(&record) {
            true => "rskip",
            false => "rkeep",
        };
        let isrev = record.is_reverse();
        let rname = std::str::from_utf8(record.qname()).unwrap();
        let mm = get_mm("MM", &record);
        let mm = match mm {
            "" => get_mm("Mm", &record),
            _ => mm,
        };
        /* mods has all the methods of an iterator,
           with item = Result<(i32 (pos), hts_base_mod)>
           plus extra methods to extract info
           about how many and what type of
           modifications are recorded fot the record
           I am looking at
        */
        let mods = match record.basemods_iter() {
            Ok(mods) => mods,
            Err(_) => {
                println!("Failed to parse MM/ML aux tags in {}\n", rname);
                continue;
            }
        };
        let mcodes: &[i32] = mods.recorded();
        let mut implicit: bool = true;
        let mut canbase: char = 'X';
        // record value of implicit for m
        for i in 0..mcodes.len() {
            let modtype = mods.query_type(mcodes[i]).unwrap();
            canbase = modtype.canonical as char;
            if (canbase == 'C' || canbase == 'G')
                && (mcodes[i] as u8 as char == 'm')
                && (modtype.implicit == 0)
            {
                implicit = false;
                break;
            };
        }
        // here depending on implicit I can iterate in different ways
        let mut calls = vec![];
        if implicit {
            call_explicit(mods, &mut calls, isrev);
            fill_implicit(&record, &mut calls, canbase);
        } else {
            call_explicit(mods, &mut calls, isrev);
        }
        sort_print_calls(mm, rskip, isrev, mcodes, implicit, rname, &mut calls);
    }
}

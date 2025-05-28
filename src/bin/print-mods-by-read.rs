use rust_htslib::htslib::BAM_FDUP;
use rust_htslib::htslib::BAM_FQCFAIL;
use rust_htslib::htslib::BAM_FSECONDARY;
use rust_htslib::htslib::BAM_FSUPPLEMENTARY;
use rust_htslib::htslib::BAM_FUNMAP;

use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;

use asms::matrix;

fn main() {
    let bamfn = "/home/ecmra/data/scan-chr6/AY7669-chr15.bam";
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    let contig = "chr15";
    let start = 1_000_000;
    let end   = 50_000_000;
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
        let ap = record.aligned_pairs();
	    let mods = match record.basemods_iter() {
		Ok(mods) => mods,
		Err(_) => {
		    println!("Failed to parse MM/ML aux tags in {}\n", rname);
		    continue
		}
	    };
	for res in mods {
            if let Ok((position, m)) = res {
		println!("{rname}:{position}:{0}:{1}:{2}",
			 m.modified_base as u8 as char,
			 m.qual, skip);
	    }
	}
    }
}

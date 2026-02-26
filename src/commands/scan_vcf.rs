use crate::utils;

use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::tbx;

use std::time;

pub fn scan_vcf(args: &[String]) {
    let usage = "usage: scan-vcf [OPTIONS] <vcf> <methfn>       
<vcf>            VCF/BCF file containing variants to look at
<methfn>         compressed (.gz) tabix indexed methylbed file
OPTIONS:
--lower-bound    <float in [0,1] default=0.4>
--upper-bound    <float in [0,1] default=0.6>
                 intermediate regions are defined as having methylation
                 between lower bound and upper bound.
--region         <CONTIG:NUM-NUM>
                 only look at specific region
";
    if (0 == args.len()) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    let mut lb:f32 = 0.4;
    let mut ub:f32 = 0.6;
    let mut region: &str = "";
    let mut i: usize = 0;
    loop {
        if i >= args.len() {
            eprintln!("{}", usage);
            return;
        }
        if "--region" == args[i] {
            region = &args[i + 1];
            i = i + 2;
            continue;
        }
	if "--lower-bound" == args[i] {
            lb = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
	if "--upper-bound" == args[i] {
            ub = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
        break;
    }
    if args.len() - i != 2 {
        eprintln!("{}", usage);
        return;
    }
    let now = time::Instant::now();
    let (vcf, methfn) = (&args[i], &args[i + 1]);
    eprintln!("vcf:{}", vcf);
    eprintln!("methfn:{}", methfn);
    let mut tbxreader =
        tbx::Reader::from_path(&methfn).expect(&format!("Could not open {}", methfn));
    let mut vcfreader = IndexedReader::from_path(vcf).expect("Error opening vcf file");
    let header_view = vcfreader.header().clone();
    if region != "" {
        let (contig, start, end) = utils::parse_region(region);
        let rid = header_view.name2rid(contig.as_bytes()).unwrap();
        vcfreader
            .fetch(rid, start as u64, Some(end as u64))
            .unwrap();
    }
    for record in vcfreader.records() {
        let record = record.unwrap();
	let lr = 2 == record.alleles().len();
        let l1 = 1 == record.alleles()[0].len();
        let l2 = 1 == record.alleles()[1].len();
	let gt = record.genotypes().expect("error in reading GT").get(0);
	let l3 = (&format!("{gt}") != "1/1") && (&format!("{gt}") != "1|1");
	let l4 = record.qual() >= 20.0; 
	if lr && l1 && l2 && l3 && l4 {
	    let contig = String::from_utf8(
                header_view
                    .rid2name(record.rid().unwrap())
                    .unwrap()
                    .to_vec(),
            )
		.unwrap();
	    let a1 = record.alleles()[0][0] as char ;
	    let a2 = record.alleles()[1][0] as char ;
            let pos = record.pos() as u64 + 1;
            let msmaybe = utils::get_bedmethyl_slice(&mut tbxreader, &contig, pos - 100, pos + 100);
	    let ms = match msmaybe {
		Ok(ms) => ms,
		Err(_) => continue
	    };
            let meanmeth: f32 = ms.meth.iter().sum::<f32>() / ms.meth.len() as f32;
            if ms.depth.len() >= 3 && meanmeth >= lb && meanmeth <= ub {
		println!(
                    "{contig}\t{}\t{}\t{pos}\t{a1}\t{a2}\t{:.4}\t{}",
                    pos - 100,
                    pos + 100,
                    meanmeth,
                    ms.depth.len(),
                );
            }
        }
    }
    let elapsed = now.elapsed();
    eprintln!("[scan-vcf] elapsed:{}sec(s)", elapsed.as_secs_f32());
}

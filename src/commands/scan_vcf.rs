use crate::utils;

use rust_htslib::bcf::{IndexedReader, Read, Reader};
use rust_htslib::tbx;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::errors::Result;
use std::time;

pub fn scan_vcf(args: &[String]) {
    let usage = "usage: scan-vcf [OPTIONS] <vcf> <methfn>       
<vcf>            VCF/BCF file containing variants to look at
<methfn>         compressed (.gz) tabix indexed methylbed file
OPTIONS:
--lower-bound    <float in [0,1] default=0.4>
--upper-bound    <float in [0,1] default=0.6>
                 intermediate regions are defined as having methylation
                 between lower bound and upper bound

--region         <CONTIG:NUM-NUM>
                 only look at specific region

--whl            half size of the window centered around the SNP
                 which is used to compute average methylation
                 and check whether it belongs to [lower bound, upper bound]
";
    if (0 == args.len()) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    let mut lb:f32 = 0.4;
    let mut ub:f32 = 0.6;
    // window half size
    let mut whl:usize = 100;
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
	if "--whl" == args[i] {
            whl = args[i + 1].parse().unwrap();
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
    
    match (vcf.as_str(), region) {
	("-", "") => {
	    let mut rdr = Reader::from_stdin().expect("could not read vcf from stdin");
	    let hv = rdr.header().clone();
	    process_vcf(&mut rdr, &hv, &mut tbxreader, lb, ub, whl)
	}
	("-", _) => panic!("can't specify region when reading from stdin"),
	(p, "") => {
	    let mut rdr = Reader::from_path(p).expect(&format!("could not open {}", p));
	    let hv = rdr.header().clone();
	    process_vcf(&mut rdr, &hv, &mut tbxreader, lb, ub, whl)
	}
	(p, r) => {
	    let mut rdr = IndexedReader::from_path(p).expect(&format!("could not open {}", p));
	    let hv = rdr.header().clone();
	    let (contig, start, end) = utils::parse_region(r);
	    let rid = hv.name2rid(contig.as_bytes()).expect(&format!("contig {} not in vcf header", contig));
	    rdr.fetch(rid, start as u64, Some(end as u64)).expect(&format!("could not fetch region {}", r));
	    process_vcf(&mut rdr, &hv, &mut tbxreader, lb, ub, whl)
	}
    }
    .expect("error while processing vcf");
    let elapsed = now.elapsed();
    eprintln!("[scan-vcf] elapsed:{}sec(s)", elapsed.as_secs_f32());
}


fn process_vcf<R:Read>(vcf: &mut R, header_view:&HeaderView, tbxreader: &mut tbx::Reader, lb:f32, ub:f32, whl:usize) ->
    Result<(), Box<dyn std::error::Error>> {
    for record in vcf.records() {
	let record = record?;
	let lr = 2 == record.alleles().len();
        let l1 = lr && 1 == record.alleles()[0].len();
        let l2 = lr && 1 == record.alleles()[1].len();
	let gt = match record.genotypes() {
	    Ok(gts) => gts.get(0),
	    Err(_) => continue
	};
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
	    let wstart = pos.saturating_sub(whl as u64);
	    let wend = pos + whl as u64;
            let msmaybe = utils::get_bedmethyl_slice(tbxreader,
					   &contig,
					   wstart,
					   wend);
	    let ms = match msmaybe {
		Ok(ms) => ms,
		Err(_) => continue
	    };
            let meanmeth: f32 = ms.meth.iter().sum::<f32>() / ms.meth.len() as f32;
            if ms.depth.len() >= 3 && meanmeth >= lb && meanmeth <= ub {
		println!(
                    "{contig}\t{}\t{}\t{pos}\t{a1}\t{a2}\t{:.4}\t{}",
                    wstart,
                    wend,
                    meanmeth,
                    ms.depth.len(),
                );
            }
        }
    }
	Ok(())
    }
	


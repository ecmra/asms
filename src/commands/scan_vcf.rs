use crate::utils;

use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::tbx;

/*
../scan-chr6/out/vcfs/tmp/phase_output/phase_vcf/phased_chr15.vcf.gz
*/

pub fn scan_vcf(args: &[String]) {
    let usage = "usage: scan-vcf [OPTIONS] <vcf> <methfn> <outdir>         
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
    let (vcf, methfn) = (&args[i], &args[i + 1]);
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
        if lr && l1 && l2 {
            let contig = String::from_utf8(
                header_view
                    .rid2name(record.rid().unwrap())
                    .unwrap()
                    .to_vec(),
            )
            .unwrap();
            let pos = record.pos() as u64 + 1;
            let ms = utils::get_bedmethyl_slice(&mut tbxreader, &contig, pos - 100, pos + 100);
            let meanmeth: f32 = ms.meth.iter().sum::<f32>() / ms.meth.len() as f32;
            if ms.depth.len() >= 3 && meanmeth >= lb && meanmeth <= ub {
                println!(
                    "{}\t{}\t{}\t{:.4}\t{}\t{}",
                    contig,
                    pos - 100,
                    pos + 100,
                    meanmeth,
                    ms.depth.len(),
                    pos
                );
            }
        }
    }
}

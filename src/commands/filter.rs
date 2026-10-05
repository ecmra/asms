//! use statistics on the clusters to call asm regions

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;

//use crate::commands::cluster;
use crate::commands::pval;
//use fxhash::FxHashMap;
//use std::path::Path;

use crate::utils;

/// read summary stats and use them to
/// test the significance of the clustering
pub fn filter(args: &[String]) {
    let usage = "usage: filter <outdir>
<outdir>         directory which contains the stats subdirectory.
                 the output of filter is in <outdir>/pass.bed 
";
    if ( args.len() < 1 ) || ( args[0] == "-h" ) || ( args[0] == "--help" ) {
        eprintln!("{}", usage);
        return;
    }
    let outdir = &args[0];
    let asmfn = &format!("{}/pass.bed", outdir);
    let handleasm = BufWriter::new(File::create(asmfn).expect(&format!("==can't create {}==\n", asmfn)));
    print_asm(outdir, handleasm);
}

/// read summary stats and call asms
///
/// see BED standard https://samtools.github.io/hts-specs/BEDv1.pdf
/// The positions in BED
/// ﬁelds are all described in the 0-based,
/// half-open coordinate system
fn print_asm<W: std::io::Write>(outdir: &str, mut handle: BufWriter<W>) {
    let summaryfn = &format!("{}/stats-summary.bed", outdir);
    let bufreader = BufReader::new(File::open(summaryfn).expect("==can't open summary==\n"));
    utils::write_out(
        &mut handle,
        &format!("#chrom\tchromStart\tchromEnd\tname\tscore(nCpG)\tpop0\tpop1\tm0\tm1\tr\tn\tabsmdiff\tstatsfn\n"),
    );
    let pvalthr = 1e-3;
    for line in bufreader.lines() {
        let line = line.unwrap();
        if line[0..1] == *"#" {
            continue;
        }
        let fields: Vec<&str> = line.trim_end().split_whitespace().collect::<Vec<_>>();
        let chrom = &fields[0];
        let chromstart = &fields[1];
        let chromend = &fields[2];
        let pop0: usize = fields[6]
            .parse()
            .expect(&format!("can't parse {}", fields[6]));
        let pop1: usize = fields[7].parse().unwrap();
        let n: usize = fields[3].parse().unwrap();
        let d: usize = fields[5].parse().unwrap();
        let absmdiff: f32 = fields[8].parse().unwrap();
	let statsfn   = fields[11];
	let clusterfn = &format!("{}/{}", outdir, fields[12]);
	let matrixfn  = &format!("{}/{}", outdir, fields[13]);
	println!("{}", clusterfn);
	let (_refmedian, pval) = pval::pval_of_filenames(matrixfn, clusterfn);
        let r = pop0 as f32 / (pop0 + pop1) as f32;
        let c1 = ( r >= 0.2 ) && ( r <= 0.8 );
        let c2 = ( n >= 10 ) && ( d >= 2 ) && ( absmdiff >= 0.2 );
	let status:&str;
	if c1 && c2 && (pval >= 0.0 ) && ( pval <= pvalthr) {
	    status = "PASS";
	} else {
	    status = "FAIL"
	}
        utils::write_out(
                &mut handle,
                &format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{pval}\t{}\n",
                    chrom, chromstart, chromend, status, d, pop0, pop1,
		    fields[9], fields[10],r,n,absmdiff, statsfn
                ),
            );
    }
}

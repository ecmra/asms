extern crate asms;
use asms::commands;
use std::env;

const BUILD_TIME : &str = include!(concat!(env!("OUT_DIR"), "/timestamp.txt"));
const VERSION: &str = env!("CARGO_PKG_VERSION");


fn main() {
    eprintln!("== asms ==");
    eprintln!("@github.com:ecmra/asms.git");
    eprintln!("version:{}", VERSION);
    let commit = option_env!("COMMIT").unwrap_or("not defined");
    eprintln!("This binary was built from commit:[{}]", commit);
    eprintln!("compiled on {}", BUILD_TIME);
    let usage =   "\nusage: asms COMMAND
where COMMAND can be:

scan         find intermediate methylation regions in a bedmethylfile
cluster      cluster reads according to methylation patterns
filter       list putative allele specific methylation loci
scan-vcf     check for intermediate methylation around genomic variants
pval         compute p-value of clustering according to a permutation test
";
    let args:Vec<String> = env::args().collect();
    if (1 == args.len()) || (args[1] == "-h") || (args[1] == "--help")  {
	eprintln!("{}", usage);
	return;
    }
    if "scan" == args[1] {
	commands::scan::scan(&args.get(2..).unwrap());
	return
    }
    if "cluster" == args[1] {
	commands::cluster::cluster(&args.get(2..).unwrap());
	return
    }
    if "filter" == args[1] {
	commands::filter::filter(&args.get(2..).unwrap());
	return
    }
    if "scan-vcf" == args[1] {
	commands::scan_vcf::scan_vcf(&args.get(2..).unwrap());
	return
    }
    if "pval" == args[1] {
	commands::pval::pval(&args.get(2..).unwrap());
	return
    }
    eprintln!("don't recognize:{}", args[1]);
    eprintln!("{}", usage);
}

    

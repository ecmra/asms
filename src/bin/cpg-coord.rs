
use std::env;
use rust_htslib::faidx;

fn main() {
    let args:Vec<String> = env::args().collect();
    let reffn = &args[1];
    let refreader = faidx::Reader::from_path(reffn).unwrap();
    let nseq = refreader.n_seqs();
    eprintln!("{nseq} seqs found in {reffn}");
    for i in 0..nseq {
	let seqname = refreader.seq_name(i as i32).unwrap();
	let seql = refreader.fetch_seq_len(&seqname);
	eprintln!("{seqname}\t{seql}");
    }
}


use std::io::BufRead;
use std::io::BufReader;
use std::fs::File;

use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read;
use rust_htslib::htslib::{BAM_FDUP,BAM_FQCFAIL,
			  BAM_FSECONDARY,BAM_FSUPPLEMENTARY,BAM_FUNMAP};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Seq;

use fastrand;

use crate::parsemod;

pub fn cluster_by_snp(args: &[String]) {
    let usage = "usage: cluster-by-snp [OPTIONS] <bamfn> <varfn>
<bamfn>          alignment file
<varfn>          list of het snps to check

OPTIONS:
--seed           <u64>
                 seed for the random number generator 
                 (to initialize permutation testing)

--nperm          <usize>
                 number of permutations (default=1000)

";
    if (args.len() < 2) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    let mut seed:u64 = 0;
    let mut nperm:usize = 1000;
    let mut i:usize = 0;
    loop {
        if i >= args.len() {
            eprintln!("{}", usage);
            return;
        }
        if "--seed" == args[i] {
            seed = args[i + 1].parse().expect("can't parse seed");
            fastrand::seed(seed);
            i = i + 2;
            continue;
        }
        if "--nperm" == args[i] {
            nperm = args[i + 1].parse().expect("can't parse maxit");
            i = i + 2;
            continue;
        }
	break;
    }
    if args.len() - i != 2 {
        eprintln!("{}", usage);
        return;
    }
    let bamfn = &args[i];
    let varfn = &args[i + 1];
    let mut bamreader = IndexedReader::from_path(bamfn).
	expect("can't open the BAM file");
    let stdin = std::io::stdin();
    let varf:Box<dyn BufRead>;
    if "-" == varfn {
    	varf = Box::new(stdin.lock());
    } else {
    	varf = Box::new(BufReader::new(File::open(varfn).unwrap()));
    }
    let mut lc = 0;
    for line in varf.lines() {
    	let line = line.unwrap();
	lc += 1;
    	let fields:Vec<&str> = line.split(':').collect();
	let contig:&str = fields[0];
	let gpos1snp:i64 = fields[1].parse::<i64>().expect("invalid snp pos");
	let refnuc:u8 = fields[2].bytes().nth(0).expect("invalid refnuc");
	let altnuc:u8 = fields[3].bytes().nth(0).expect("invalid altnuc");
	let ( cl, mmeth, ncpgs, m01, n01 ) =
	    meth_by_snp(&mut bamreader, contig, gpos1snp, refnuc, altnuc);
	let (stat, pval) = test_snp_perm(cl.into(), mmeth, ncpgs, nperm, m01, n01);
	println!("{lc}:{line}:{pval:.3e}:{:.2}:{}:{:.2}:{}", m01[0],n01[0],m01[1],n01[1]);
    }
}

fn meth_by_snp(bamreader:&mut IndexedReader, contig:&str, gpos1snp:i64, refnuc:u8, altnuc:u8) ->
    (Vec<usize>, Vec<f32>, Vec<usize>, [f32;2], [usize;2]) {
        let mut recordbuf = Record::new();
	// methylation along the read
	let mut meth:Vec<f32>;
	let start:i64 = gpos1snp - 500 ;
	let end:i64 = gpos1snp + 500 ;
	bamreader.fetch((contig, start, end)).unwrap();
	let mut cl:Vec<usize> = Vec::with_capacity(1000); 
	let mut mmeth:Vec<f32> = Vec::with_capacity(1000);
	let mut ncpgs:Vec<usize> = Vec::with_capacity(1000);
	let mut rc:usize=0; // read counter
	let mut m01:[f32;2] = [0.,0.];
	let mut n01:[usize;2] = [0,0];
	while let Some(result) = bamreader.read(&mut recordbuf) {
	    let record = match result {
		Ok(_) => &recordbuf,
		Err(_) => panic!("BAM parsing failed"),
	    };
	    let flags = record.flags();
	    let qflags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY
		| BAM_FQCFAIL | BAM_FDUP; 
	    if flags & qflags as u16 != 0 { continue }
	    if record.mapq() < 10 { continue }
	    let rname = std::str::from_utf8(record.qname()).unwrap();
	    let seq = record.seq();
	    let ap =  record.aligned_pairs();
	    let (mut minrpos0, mut maxrpos0, mut rpos0snp): (i64, i64, usize) =
		(i64::MAX, 0, 0);
	    let mut found:bool = false;
	    for [rpos0, gpos0] in ap {
		let gpos1 = gpos0+1;
		if gpos1 < start || gpos1 > end {continue}
		if rpos0 < minrpos0 {minrpos0=rpos0}
		if rpos0 > maxrpos0 {maxrpos0=rpos0}
		if gpos1 == gpos1snp {
		    rpos0snp = rpos0 as usize; found = true;
		}
	    }
	    if !found {continue} 
	    let nuc = unsafe{ seq.decoded_base_unchecked(rpos0snp)};
	    let hp = match nuc {
		_ if nuc == refnuc => 0usize,
		_ if nuc == altnuc => 1usize,
		_ => continue,
	    };
	    meth = vec![];
	    let mut calls = vec![];
	    parsemod::parse_mods(&record, &mut calls, seq, minrpos0 as i32, maxrpos0 as i32);
	    calls.sort_by(|c1, c2| c1.0.cmp(&c2.0));
	    // calls contain meth info for all the Cs not specifically the CGs.
	    for res in calls {
		let modrpos0 = res.0 as usize;
		let m = res.1;
		let (nuc0, nuc1) = match record.is_reverse() {
		    false => unsafe {base_pair_at_rpos0(&seq, modrpos0)},
		    true => unsafe{base_pair_at_rpos0(&seq, modrpos0-1)} 
		};
		// C and G could be nearby on the read, and far away on the genome
		// because of a deletion; I assume this does not happen frequently,
		// hence I don't check.
		let checkcg = ( nuc0 == b'C' ) && ( nuc1 == b'G' );
		if !(b'm' == m.modified_base as u8) || !checkcg {continue}
		let mqual = m.qual as u8;
		if mqual > 127 {meth.push(1.0)}
		else {meth.push(0.0)}
	    }
	    let ncg = meth.len();
	    let mm = meth.iter().sum::<f32>()/ncg as f32;
	    let sign = match record.is_reverse() {
		false => '+',
		true  => '-',
	    };
	    let (rn, an) = (refnuc as char, altnuc as char);
	    /*match hp {
		0 => println!("{rname}\t0\t{mm}\t{ncg}\t{rn}>{an}\t{sign}"),
		1 => println!("{rname}\t1\t{mm}\t{ncg}\t{rn}>{an}\t{sign}"),
		_ => panic!("EEEEEEEEEEEEE????"),
	    }*/
	    cl.push(hp); mmeth.push(mm); ncpgs.push(ncg);
	    n01[cl[rc]] +=1;
	    m01[cl[rc]] += mmeth[rc];
	    rc += 1;
	} // while loop on reads 
	m01[0] = m01[0] / n01[0] as f32;
	m01[1] = m01[1] / n01[1] as f32;
    	return (cl, mmeth, ncpgs, m01, n01);
    }

unsafe fn base_pair_at_rpos0(seq:&Seq<'_>, rpos0:usize) -> (u8,u8) {
    return (unsafe{ (*seq).decoded_base_unchecked(rpos0) }, 
     unsafe{ (*seq).decoded_base_unchecked(rpos0 + 1) })
}

/// test if abs(m0 - m1) is significantly different from 0
///
/// the background distribution is built by
/// randomizing the cluster assignment
fn test_snp_perm(cl:Vec<usize>, mmeth:Vec<f32>, ncpgs:Vec<usize>,
		 nperm:usize, m01:[f32;2], n01:[usize;2])
	    -> (f32, f32) {
    let refdiff = f32::abs(m01[0] - m01[1]);
    let mut m01:[f32;2] = [0.,0.];
    let mut n01:[usize;2] = [0,0];
    let mut ngreater:usize = 0;
    for p in 0..nperm {
	for i in 0..cl.len() {
	    let idx = fastrand::usize(0..2);
	    n01[idx]+=1; m01[idx] += mmeth[i];
	}
	m01[0] = m01[0] / n01[0] as f32;
	m01[1] = m01[1] / n01[1] as f32;
	let rnddiff = f32::abs(m01[0] - m01[1]);
	if rnddiff > refdiff {ngreater+=1}
    }
    return (ngreater as f32 , ngreater as f32 / nperm as f32);
    //    return Err("oh my!".to_string().into());
}

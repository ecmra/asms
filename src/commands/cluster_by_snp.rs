//! cluster by read can print the cluster assignment for each read.
//! 
//! - it takes in input (optionally) a file name
//!   where to write the cluster assignment for each input read
//! - if no file name is given, the program behaves as before,
//!   and prints no clustering information. Otherwise, see below
//! - the cluster assignment is 0 for ref and 1 for alternative
//! - all reads are printed in the same file, which has 
//!   the following output 4 columns:
//!   1 rname
//!   2 cluster 0/1
//!   3 snp coordinates+ref+alt
//!   4 avg methylation in window
//! - at the bottom of the file there's a set
//!   of summary statistics, prefixed by #@
//!   #@N total number of reads
//!   #@WS window size to compute average methylation
//!   #@POP0/POP1 number of ref/alt reads
//!   #@DT date/time



use std::io::BufRead;
use std::io::BufReader;
use std::fs::File;
use std::time;

use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read;
use rust_htslib::htslib::{BAM_FDUP,BAM_FQCFAIL,
			  BAM_FSECONDARY,BAM_FSUPPLEMENTARY,BAM_FUNMAP};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Seq;

use fastrand;

use crate::parsemod;

use numutils;


// struct Snp<'a>{
//     contig: &'a str,
//     gpos1snp: i64,
//     refnuc: u8,
//     altnuc: u8,
// }

/// filter vectors positions looking at the corresponding index in ncpgs
///
/// use to filter away reads which do not contain enough CpGs
/// in the measurement window
fn filter_by_ncpgs<T:Copy>(v:&[T], ncpgs:&[usize], thr:usize) -> Vec<T>   {
    return (*v).iter().
    enumerate().
    filter(|x| ncpgs[(*x).0]> thr).
    map(|x| *(x.1)).collect();
}

/// compute mean methylation by cluster
///
/// fill two array:
///    m01 contains the means
///    n01 contains the number of elements per cluster
pub fn mmeth_by_cluster(cl:&[usize], mmeth:&[f32]) -> ([f32; 2], [usize; 2]) {
    let mut  m01:[f32;2]= [0.0,0.0]; 
    let mut n01:[usize;2] = [0,0];
    for i in 0..cl.len() {
	let idx = cl[i];
	m01[idx]+=mmeth[i]; n01[idx]+=1;
    }
    m01[0]=m01[0]/n01[0] as f32;
    m01[1]=m01[1]/n01[1] as f32;
    return (m01, n01);
}

pub fn cluster_by_snp(args: &[String]) {
    let usage = "usage: cluster-by-snp [OPTIONS] <bamfn> <varfn>
<bamfn>          alignment file
<varfn>          list of het snps to check

OPTIONS:
--seed           <u64>
                 seed for the random number generator 
                 (to initialize permutation testing)

--maxperm          <usize>
                 max number of permutations (default=1000)

";
    if (args.len() < 2) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    let mut seed:u64;
    let mut maxperm:usize = 1000;
    let mut i:usize = 0;
    //let mut clusterfn:Option<&str> = None;
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
            maxperm = args[i + 1].parse().expect("can't parse maxperm");
            i = i + 2;
            continue;
        }
	/*if "--clusterfn" == args[i] {
            clusterfn = Some(&args[i + 1]);
            i = i + 2;
            continue;
        }*/
	break;
    }
    if args.len() - i != 2 {
        eprintln!("{}", usage);
        return;
    }
    let now = time::Instant::now();
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
    let mut lines:Vec<String> = vec![];
    let mut pvals:Vec<f32> = vec![];
    for line in varf.lines() {
	let line = line.unwrap();
	eprintln!("{}", line);
	lc = lc + 1;
    	let fields:Vec<&str> = line.split('\t').collect();
	let contig:&str = fields[0];
	let gpos1snp:i64 = fields[3].parse::<i64>().expect("invalid snp pos");
	let refnuc:u8 = fields[4].bytes().nth(0).expect("invalid refnuc");
	let altnuc:u8 = fields[5].bytes().nth(0).expect("invalid altnuc");
	let _crefnuc = refnuc as char;
	let _caltnuc = altnuc as char;
	let ( cl, mmeth, ncpgs, _m01, _n01 ) =
	    meth_by_snp(&mut bamreader, contig, gpos1snp, refnuc, altnuc);
	let cl:Vec<usize> = filter_by_ncpgs(&cl, &ncpgs, 2); 
        let mmeth:Vec<f32> =  filter_by_ncpgs(&mmeth, &ncpgs, 2); 
        let _ncpgs:Vec<usize> = filter_by_ncpgs(&ncpgs, &ncpgs, 2);
	let (m01, n01) = mmeth_by_cluster(&cl, &mmeth); 
	let (_stat, pval, refdiff, _nperm ) =
	    test_snp_perm(&cl, &mmeth, maxperm, m01, n01, 1u8);
	// meth in cluster 0, #reads in cluster 0, meth in cluster 1, #reads in cluster 1
	let rec = format!("{line}\t{pval:.3e}\t{:.2}\t{}\
			   \t{:.2}\t{}\t{:.2}",
			  m01[0],n01[0],m01[1],n01[1], refdiff);
	lines.push(rec); pvals.push(pval);
    }
    print_with_adj_pval(&lines, &pvals);
    let elapsed = now.elapsed();
    eprintln!("[cluster-by-snp] elapsed:{}sec(s)", elapsed.as_secs_f32());
}

fn print_with_adj_pval(lines:&Vec<String>, pvals:&Vec<f32>){
    let adjpvals = numutils::ch14::padjust(pvals);
    let nlines:usize = lines.len();
    for i in 0..nlines {
	println!("{}\t{:.3e}", lines[i], adjpvals[i]);
    }
}




/// cluster reads based on snp
///
/// returns:
/// Vec<usize> cluster for read i where 0=ref, 1=alt
/// Vec<f32> mean methylation for read i in window
/// Vec<usize> n cpgs on read i in wondows
/// [f32;2] mean methylation in cluster 0/1
/// [usize;2] number of reads in cluster 0/1
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
	    let _rname = std::str::from_utf8(record.qname()).unwrap();
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
		// when the read is reversed the coordinate on the
		// read modrpos0 corresponds to a G.
		// when checking that I am looking at a CG dinucleotide
		// and the read is reversed I substract -1 from modrpos0 so that
		// it points to the adjacent C.
		// I can't do that when the read is revered and modrpos0 == 0.
		// If this is the case, I skip the call.
		if (modrpos0==0) && (record.is_reverse()) {continue;} 
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
	    let _sign = match record.is_reverse() {
		false => '+',
		true  => '-',
	    };
	    let (_rn, _an) = (refnuc as char, altnuc as char);
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


/// return nucleotide pair at given position in the read
unsafe fn base_pair_at_rpos0(seq:&Seq<'_>, rpos0:usize) -> (u8,u8) {
    return (unsafe{ (*seq).decoded_base_unchecked(rpos0) }, 
     unsafe{ (*seq).decoded_base_unchecked(rpos0 + 1) })
}


/// test if abs(m0 - m1) is significantly different from 0
///
/// the background distribution is built by
/// randomizing the cluster assignment
/// == arguments:
/// mmeth[i]: mean methylation in window for read i
/// ncpgs[i]: number of cpgs in window for read i
/// maxperm: number of permutations
/// m01: mean methylation in cluster 0/1
/// n01:number of reads in cluster 0/1
/// mode :u8  0  for vanilla Montecarlo, 1 for sequential Montecarlo
///        the reference for sequential Montecarlo is
///        Besag/Clifford Sequential montecarlo p values (Biometrika 1991)
/// == return:
/// ngreater: usize number of random values greater than the reference
///           (real) value
/// pvalue:   f32
/// refdiff:  f32
/// nperm:    usize number of permutations acutally sampled
///           (can be lower than maxperm in sequential mode )
pub fn test_snp_perm(cl:&Vec<usize>, mmeth:&Vec<f32>, 
		 maxperm:usize, m01:[f32;2], n01:[usize;2],mode:u8)
	    -> (usize, f32, f32, usize) {
    let refdiff = f32::abs(m01[0] - m01[1]);
    let mut m01:[f32;2];
    let mut ngreater:usize = 0;
    let mut cl = cl.clone();
    // if early stopping I stop after I see
    // h values greater than refdiff
    let h = 20; 
    // number of permutations sampled
    let mut nperm:usize=0;
    for _ in 0..maxperm {
	// number of reads and methylation for allele 0,1
	m01 = [0.,0.];
	// shuffle preserve the cardinality of the clusters
	fastrand::shuffle(&mut cl);
	nperm += 1;
	for i in 0..cl.len() {
	    let idx = cl[i];
	    m01[idx] += mmeth[i];
	}
	m01[0] = m01[0] / (n01[0] as f32);
	m01[1] = m01[1] / (n01[1] as f32);
	let rnddiff = f32::abs(m01[0] - m01[1]);
	if rnddiff > refdiff { ngreater += 1 }
	if ngreater == h && (mode == 1)  {
	    return (ngreater  , ngreater as f32 / nperm  as f32, refdiff, nperm);
	} 
    }
    // this is the same expression regardless of whether the function
    // is called in mode 0 or mode 1
    // mode 0: always reaches here
    // mode 1: reaches here only if ngreater never hit h (low p-value regime)
    return (ngreater  , (ngreater + 1) as f32 / (nperm + 1) as f32, refdiff, nperm);
}


#[test]
fn test_cluster_by_snp(){
    let cl = vec![0,0,0,1,1,1];
    let mmeth = vec![0.1,0.1,0.1,0.7,0.7,0.7];
    let maxperm = 1000;
    let m01 = [0.1, 0.7];
    let n01 = [3,3];
    let (ngreater, _prop, refdiff, _nperm) = test_snp_perm(&cl, &mmeth,  maxperm, m01, n01,0);
    assert!(ngreater==0);
    assert!(f32::abs(refdiff - 0.6) < 1e-6);
}

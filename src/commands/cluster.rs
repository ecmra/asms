//! clusters reads mapping over a predefined set of loci looking
//! at their methylation patterns.


use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::time;


use rust_htslib::bam::IndexedReader;


use rust_htslib::faidx;


use libc;



use fastrand;


use crate::utils;
use crate::stats;
use crate::matrix;
use crate::prob;


/// stores coordinates of genomic region
pub struct Region<'a> {
    pub contig: &'a str,
    pub start: i64,
    pub end: i64,
}


/// cluster data
///
/// cprobs prob for a read to belong to cluster 0/1 (sum to 1)
/// clusters for each read the cluster which has higher prob
/// mu (k X d matrix, stored in rows) mean of each cluster across positions
/// pop number of reads assigned to cluster 0/1
pub struct Cluster {
    pub cprobs: Vec<[f32; 2]>,
    pub clusters: Vec<u8>,
    pub mu: Vec<f32>,
    pub pop: [f32; 2],
}

/// clusters reads in the loci specified in `bedfn`.
///
/// Alignment and methylation information is taken from `bamfn`.
/// `reference` is used to check for the presence of CG dinnucleotides
/// in the reads and results are write in subdirectories of `outdir`.
/// possible to cluster on a set of precomputed matrices)
/// that is without recomputing the matrices.
/// cluster --precomputed-matrix  bedfile outdir looks in outdir/matrices
/// for the input to the clustering algorithms
/// checking against the bedfile to compute the names of the output files.
pub fn cluster(args: &[String]) {
    let usage = "usage: cluster [OPTIONS] <bedfn> <bamfn> <reference> <outdir>
<bedfn>          BED files containing regions for the clustering algorithm to work on
<bamfn>          alignment file
<reference>      indexed FASTA reference
<outdir>         directory which will contain the results

OPTIONS:
--seed           <u64>
                 seed for the random number generator 
                 (to initialize clustering)

--maxit          <u16>
                 max number of EM iterations (default=10)

";
    if (args.len() < 4) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    let mut seed: u64;
    let mut i:usize = 0;
    let mut maxit:u16 = 10;
    loop {
        if i >= args.len() {
            eprintln!("{}", usage);
            return;
        }
        if "--seed" == args[i] {
            seed = args[i + 1].parse().unwrap();
            fastrand::seed(seed);
            i = i + 2;
            continue;
        }
        if "--maxit" == args[i] {
            maxit = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
	break;
    }
    if args.len() - i != 4 {
        eprintln!("{}", usage);
        return;
    }
    let (bedfn, bamfn, reffn, outdir) = (&args[i], &args[i + 1], &args[i + 2], &args[i + 3]);
    eprintln!("bedfn:{}", bedfn);
    eprintln!("bamfn:{}", bamfn);
    eprintln!("reffn:{}", reffn);
    eprintln!("outdir:{}", outdir);
    let now = time::Instant::now();
    fs::create_dir_all(outdir).expect(&format!("can't create {}", outdir));
    fs::create_dir_all(&format!("{}/matrices", outdir)).unwrap();
    fs::create_dir_all(&format!("{}/clusters", outdir)).unwrap();
    fs::create_dir_all(&format!("{}/stats", outdir)).unwrap();
    let refreader = faidx::Reader::from_path(reffn).unwrap();
    let fstname = refreader.seq_name(0i32).unwrap();
    let bufreader = BufReader::new(File::open(bedfn).expect("===can't open BED file==="));
    let mut bamreader = IndexedReader::from_path(bamfn).unwrap();
    let _ = bamreader.set_reference(Path::new(reffn));
    let mut oldcontig = fstname.to_string();
    let mut reflen: usize = refreader.fetch_seq_len(&fstname) as usize;
    let mut refseq: &[u8] = refreader.fetch_seq(&fstname, 0, reflen - 1).unwrap();
    for line in bufreader.lines() {
        let line = line.unwrap();
        if line.starts_with("#") {
            continue;
        }
        let fields: Vec<&str> = line.trim_end().split_whitespace().collect::<Vec<_>>();
        let contig = fields[0];
        let start: i64 = fields[1].parse::<i64>().unwrap();
        let end: i64 = fields[2].parse::<i64>().unwrap();
        let region = Region {
            contig: contig,
            start: start,
            end: end,
        };
        eprintln!("processing {}:{}-{}", contig, start, end);
        let nowpartial = time::Instant::now();
        if contig != oldcontig {
            // Free up memory
            // see the discussion at
            // https://github.com/rust-bio/rust-htslib/issues/401
            unsafe { libc::free(refseq.as_ptr() as *mut std::ffi::c_void) };
            reflen = refreader.fetch_seq_len(contig) as usize;
            eprintln!("loading:{}({}bp(s))", contig, reflen);
            refseq = refreader.fetch_seq(contig, 0, reflen - 1).unwrap();
        }
        let m = matrix::meth_of_bam(&mut bamreader, region, outdir, refseq);
        eprintln!("found {} reads, {} CpGs", m.n, m.d);
        if (0 == m.n) || (0 == m.d) {
            continue;
        }

        let outfn = &format!(
            "{}/clusters/cluster-{}-{}-{}.txt",
            outdir, contig, start, end
        );
        let mut handle = BufWriter::new(File::create(outfn).unwrap());
        let cm = cluster_matrix(&m, maxit, &mut handle);
        print_clusters(&m, &cm.cprobs, &cm.clusters, &cm.mu, &cm.pop, handle);
        let outfn = &format!("{}/stats/stats-{}-{}-{}.txt", outdir, contig, start, end);
        let handle = BufWriter::new(File::create(outfn).unwrap());
        let stats = stats::compute_stats(&m, &cm.clusters);
	stats::print_stats(&stats, handle);
	let elapsedpartial = nowpartial.elapsed();
        eprintln!("done in {}sec(s)", elapsedpartial.as_secs_f32());
        oldcontig = contig.to_string();
    }
    let summaryfn = &format!("{}/stats-summary.bed", outdir);
    let handle = BufWriter::new(File::create(summaryfn).unwrap());
    stats::print_summary(bedfn, outdir, handle);
    let elapsed = now.elapsed();
    eprintln!("cluster ran in {}sec(s)", elapsed.as_secs_f32());
}

pub fn cluster_of_file(fname:&str) -> Cluster {
    let mut clusters:Vec<u8>=vec![];
    let bufreader = BufReader::new(File::open(fname).unwrap());
    for line in bufreader.lines() {
	let s = line.unwrap();
	if &s[0..2] == "#@" {continue;}
	let fields = s.split('\t').collect::<Vec<&str>>();
	clusters.push(fields[1].parse::<u8>().unwrap());
    }
    let cl = Cluster{ cprobs:vec![[0.0f32,0.0f32]],
		      clusters: clusters,
		      mu: vec![-1.0],
		      pop: [-1.0,-1.0]
    };
    return cl;
}


/// use mu and pi to compute gamma  (matrix n times 2)
///
/// see the EM chapter in PRML for what gamma
/// means in EM.
/// gamma is stored by row into a 1d vector.
/// hence position (i,j) is stored in g\[i*2+j\].
fn bernoulli_gamma(m: &matrix::Matrix, pi: &[f32; 2], mu: &Vec<f32>)
		   -> Vec<f32> {
    let n = m.n;
    let d = m.d;
    let mut g: Vec<f32> = vec![0.0; n * 2];
    for i in 0..n {
        let obs = &m.m[i * d..d * (i + 1)];
        //println!("{:?}", obs);
        let log_gi0 = f32::ln(pi[0]) + prob::logmvbern(obs, &mu[0..d]);
        let log_gi1 = f32::ln(pi[1]) + prob::logmvbern(obs, &mu[d..2 * d]);
        let logdenom = prob::logsumexp(&[log_gi0, log_gi1]);
        g[i * 2] = f32::exp(log_gi0 - logdenom);
        g[i * 2 + 1] = f32::exp(log_gi1 - logdenom);
    }
    return g;
}

/// update mu and pop
fn update_mu(m: &matrix::Matrix, gamma: &Vec<f32>) -> (Vec<f32>, [f32; 2]) {
    let (mut pop0, mut pop1) = (0.0f32, 0.0f32);
    for i in 0..m.n {
        pop0 += gamma[i * 2];
        pop1 += gamma[i * 2 + 1];
    }
    let mut res: Vec<f32> = vec![0.0; 2 * m.d];
    for j in 0..m.d {
        for i in 0..m.n {
            let meth = m.m[i * m.d + j];
            if meth != -1.0 {
                res[j] += gamma[i * 2] * meth / pop0;
                res[m.d + j] += gamma[i * 2 + 1] * meth / pop1;
            }
        }
    }
    return (res, [pop0, pop1]);
}

/// initial value for pi built using an initial clustering.
///
/// Initial clustering can be anything, including random.
/// cluster\[i\] is the cluster for read i
fn init_pi(cluster: &Vec<u8>, n: usize) -> [f32; 2] {
    let mut pi: [f32; 2] = [0.0, 0.0];
    for i in 0..n {
        pi[cluster[i] as usize] += 1.0 / n as f32;
    }
    return pi;
}

/// initial value for mu
///
/// built using any (maybe random) clustering.
/// cluster\[i\] is the cluster for read i.
/// mu is 2 times d array represented as a 1d vector.
/// mu\[i,j\] => mu \[i*d+j\] where i=0,1 j=0..d
fn init_mu(cluster: &Vec<u8>, m: &matrix::Matrix) -> Vec<f32> {
    let (n, d) = (m.n, m.d);
    let mut mu: Vec<f32> = vec![0.0; 2 * d];
    // how many reads in cluster 0,1
    let mut n01: [usize; 2] = [0, 0];
    for i in 0..n {
        let k = cluster[i] as usize;
        n01[k] += 1;
        for j in 0..d {
            let meth = m.m[i * d + j];
            // some positions are unknown (-1) and
            // are not used in the EM algorithm
            // or its initialization.
            if meth != -1.0 {
                mu[k * d + j] += m.m[i * d + j];
            }
        }
    }
    for j in 0..d {
        mu[j] = mu[j] / n01[0] as f32;
        mu[d + j] = mu[d + j] / n01[1] as f32;
    }
    return mu;
}

//log likelihood given the clusters 
#[allow(dead_code)]
fn llcluster(cl:Cluster, m:matrix::Matrix) -> f32 {
    let mut ll:f32 = 0.0;
    for i in 0..m.n {
	let clidx = cl.clusters[i] as usize;
	ll += prob::logmvbern(&m.m[i*m.d..(i*m.d+m.d)],
			      &cl.mu[clidx*m.d..clidx*m.d+m.d]);
    }
    return ll;
}

/// cluster binary matrix using a multivariate bernoulli
///
/// see PRML for a general description of the algorithm
/// and examples on gaussian mixtures. This is like
/// a mixture of multivariate bernoulli distributions.
/// returns for each read the probability that
/// it belongs to cluster 0 or cluster 1.
fn cluster_matrix<W: std::io::Write>(
    m: &matrix::Matrix,
    maxit: u16,
    mut handle: &mut BufWriter<W>,
) -> Cluster {
    let n = m.n;
    let mut cluster: Vec<u8> = vec![0; n];
    for i in 0..n {
        cluster[i] = fastrand::u8(0..2);
    }
    let mut pi = init_pi(&cluster, n);
    let mut mu = init_mu(&cluster, &m);
    // compute and print initial ll here
    utils::write_out(
        &mut handle,
        &format!("#@INITPI:{:.3}\t{:.3}\n", pi[0], pi[1]),
    );
    // mu
    utils::write_out(&mut handle, &format!("#@INITMU0:"));
    for j in 0..(m.d - 1) {
        utils::write_out(&mut handle, &format!("{:.4}\t", mu[j]));
    }
    utils::write_out(&mut handle, &format!("{:.4}\n", mu[m.d - 1]));
    utils::write_out(&mut handle, &format!("#@INITMU1:"));
    for j in 0..(m.d - 1) {
        utils::write_out(&mut handle, &format!("{:.4}\t", mu[m.d + j]));
    }
    utils::write_out(&mut handle, &format!("{:.4}\n", mu[m.d + m.d - 1]));
    let mut nit: u16 = 0;
    let mut pop: [f32; 2];
    let mut gamma: Vec<f32> = vec![0.0; 2 * n];
    let mut cprob: Vec<[f32; 2]> = vec![[0.0, 0.0]; n];
    /*  EM loop  */
    while nit < maxit {
        gamma = bernoulli_gamma(&m, &pi, &mu);
        let mupop = update_mu(&m, &gamma);
        mu = mupop.0;
        pop = mupop.1;
        pi[0] = pop[0] / n as f32;
        pi[1] = pop[1] / n as f32;
        nit += 1;
    }
    /*  *******  */
    pop = [0.0, 0.0];
    for i in 0..n {
        cprob[i] = [gamma[2 * i], gamma[2 * i + 1]];
        if cprob[i][0] > cprob[i][1] {
            cluster[i] = 0;
            pop[0] += 1.0;
        } else {
            cluster[i] = 1;
            pop[1] += 1.0;
        }
    }
    // compute and print final LL here
    return Cluster {
        cprobs: cprob,
        clusters: cluster,
        mu: mu,
        pop: pop,
    };
}

/// print clusters  on W
///
/// also prints headers #@
fn print_clusters<W: std::io::Write>(
    m: &matrix::Matrix,
    cprobs: &Vec<[f32; 2]>,
    clusters: &Vec<u8>,
    mu: &Vec<f32>,
    pop: &[f32; 2],
    mut handle: BufWriter<W>,
) {
    utils::write_out(&mut handle, &format!("#@N:{}\n#@K:2\n#@D:{}\n", m.n, m.d));
    // mu
    utils::write_out(&mut handle, &format!("#@MU0:"));
    for j in 0..(m.d - 1) {
        utils::write_out(&mut handle, &format!("{:.4}\t", mu[j]));
    }
    utils::write_out(&mut handle, &format!("{:.4}\n", mu[m.d - 1]));
    utils::write_out(&mut handle, &format!("#@MU1:"));
    for j in 0..(m.d - 1) {
        utils::write_out(&mut handle, &format!("{:.4}\t", mu[m.d + j]));
    }
    utils::write_out(&mut handle, &format!("{:.4}\n", mu[m.d + m.d - 1]));
    // pop
    utils::write_out(&mut handle, &format!("#@POP0:{}\n", pop[0]));
    utils::write_out(&mut handle, &format!("#@POP1:{}\n", pop[1]));
    // pi
    utils::write_out(&mut handle, &format!("#@PI0:{}\n", pop[0] / m.n as f32));
    utils::write_out(&mut handle, &format!("#@PI1:{}\n", pop[1] / m.n as f32));
    // clusters
    for i in 0..clusters.len() {
        utils::write_out(
            &mut handle,
            &format!(
                "{}\t{}\t{:.4}\t{:.4}\n",
                m.rnames[i], clusters[i], cprobs[i][0], cprobs[i][1]
            ),
        );
    }
}

#[test]
fn bernoulli_gamma_test() {
    let tm1: matrix::Matrix = matrix::test::make_test_matrix();
    let pi = [0.5, 0.5];
    let mu = vec![0.2, 0.8, 0.2, 0.8, 0.2, 0.8, 0.2, 0.8, 0.2, 0.8];
    let gamma = bernoulli_gamma(&tm1, &pi, &mu);
    let c0 = f32::abs(0.999024 - gamma[0]) < 1e-6;
    let c1 = f32::abs(0.00097561 - gamma[1]) < 1e-4;
    //eprintln!("{} {}", gamma[0], gamma[1]);
    assert!(c0 && c1);
}

#[test]
fn init_mu_test() {
    let tm1: matrix::Matrix = matrix::test::make_test_matrix();
    let cluster: Vec<u8> = vec![0, 1, 0, 1];
    let mu = init_mu(&cluster, &tm1);
    assert!(mu == vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
}

#[test]
fn init_pi_test() {
    let cluster: Vec<u8> = vec![0, 1, 0, 1];
    let n: usize = 4;
    let pi = init_pi(&cluster, n);
    assert!(pi == [0.5, 0.5]);
}

#[test]
fn update_mu_test() {
    let tm1: matrix::Matrix = matrix::test::make_test_matrix();
    let gamma = vec![0.9, 0.1, 0.1, 0.9, 0.9, 0.1, 0.1, 0.9];
    let mu_pop = update_mu(&tm1, &gamma);
    let mu = mu_pop.0;
    let pop = mu_pop.1;
    let c1 = pop == [2., 2.];
    let c2 = mu == vec![0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9];
    assert!(c1 && c2);
}



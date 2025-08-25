use fastrand;
use crate::stats;
use crate::matrix;
use crate::commands::cluster;

pub fn pval(args: &[String]) {
    let usage = "usage: pval <clusterfn> <matrixfn> 
<clusterfn>         cluster file
<matrixfn>          matrix file
";
    if (args.len() < 2) || (args[0] == "-h") || (args[0] == "--help") {
        eprintln!("{}", usage);
        return;
    }
    /* not used now; there might be options to
    parse later */
    let i: usize = 0; 
    let (clusterfn, matrixfn) = (&args[i], &args[i + 1]);
    let res = pval_of_filenames( &matrixfn, &clusterfn);
    println!("{}\t{}\t{}\t{}", clusterfn, matrixfn, res.0, res.1);
}

/// empirical p-value of the median absolute difference in methylation
/// between two clusters.
///
/// matrix kept fixed, this function reshuffles the cluster labels.
pub fn pval_of_matrix(m: &matrix::Matrix, cl: &Vec<u8>, nit:usize)
	    -> (f32, f32) {
    let refstats  = stats::compute_stats(m, cl);
    let refmedian = refstats.median_abs_diff_notna_0_1; 
    let mut rndmedian:f32;
    let mut pval:usize = 0;
    if cl.len() < 8 {
	eprintln!("can't compute empirical pval, too few reads");
	return (refmedian, -1.0)
    }
    let mut rndcl = cl.to_vec();
    for _ in 0..nit {
        fastrand::shuffle(&mut rndcl);
        let rndstats = stats::compute_stats(m, &rndcl);
        rndmedian = rndstats.median_abs_diff_notna_0_1;
        if  rndmedian >= refmedian { pval += 1 }
    }
    return (refmedian, pval as f32/nit as f32)
}

pub fn pval_of_filenames(fnmatrix:&str, fncluster:&str) -> (f32, f32) {
    let m = matrix::matrix_of_file(fnmatrix);
    let cl = cluster::cluster_of_file(fncluster);
    let res = pval_of_matrix(&m, &cl.clusters, 5000);
    return res;
}



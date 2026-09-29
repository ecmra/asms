
extern crate asms;
use asms::commands;
use asms::matrix;
use asms::commands::cluster;
use std::io;

fn main() -> Result<(), io::Error> {
    let fnmatrix = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/matrices/matrix-chr15-100863815-100863868.txt";
    let fncluster = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/clusters/cluster-chr15-100863815-100863868.txt";
    let mut res:(f32,f32) = (-1.0,-1.0);
    let matrix = matrix::matrix_of_file(&fnmatrix);
    let cluster = cluster::cluster_of_file(&fncluster);
    let cl = cluster.clusters;
    for _ in 0..100 {
	//res = commands::pval::pval_of_filenames(fnmatrix, fncluster);
	res = commands::pval::pval_of_matrix(&matrix, &cl, 1000);
    }
    println!("{}\t{}", res.0, res.1);
    Ok(())
}

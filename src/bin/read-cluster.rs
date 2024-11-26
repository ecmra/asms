extern crate asms;
use asms::commands::cluster;

fn main() {
    let fncluster = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/clusters/cluster-chr15-100863815-100863868.txt";
    let mut _cluster;
    for _ in 0..10{
	_cluster = cluster::cluster_of_file(&fncluster);
    }
}

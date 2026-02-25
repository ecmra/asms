use asms::commands::cluster_by_snp;

fn main() {
    let cl = vec![0,0,0,0,0,0,0,0,0,01,1,1,1,1,1,1,1,1,1];
    let mut mmeth:Vec<f32>; 
    let ncpgs = vec![3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
    let maxperm = 10000;
    let nreads = cl.len();
    for _ in 0..10000 {
	mmeth = (0..nreads).map(|_| fastrand::f32()).collect();
	let (m01, n01) = cluster_by_snp::mmeth_by_cluster(&cl, &mmeth);
	//let refdiff = f32::abs(m01[0] - m01[1]);
	let (ngreater0, pval0, refdiff, nperm0) =
	    cluster_by_snp::test_snp_perm(&cl, &mmeth,  maxperm, m01,0);
	let (ngreater1, pval1, refdiff, nperm1) =
	    cluster_by_snp::test_snp_perm(&cl, &mmeth,  maxperm, m01,1);
	print!("{}\t{}\t{}", ngreater0, pval0, nperm0);
	println!("\t{}\t{}\t{}\t{}", ngreater1, pval1, nperm1, refdiff);
    }
}

use asms::commands::cluster_by_snp;

fn main() {
    let cl = vec![0,0,0,1,1,1];
    let mut mmeth:Vec<f32>; 
    let ncpgs = vec![3,3,3,3,3,3];
    let maxperm = 1000;
    let m01 = [0.1, 0.7];
    for _ in 0..1000 {
	mmeth= (0..6).map(|_| fastrand::f32()).collect();
	for i in 0..3 {
	    mmeth[i] = mmeth[i] - 0.4;
	    if mmeth[i] <= 0.0 { mmeth[i]=0.0}
	}
	eprintln!("{:?}", mmeth);
	let (ngreater0, pval0, refdiff, nperm0) =
	    cluster_by_snp::test_snp_perm(&cl, &mmeth, &ncpgs, maxperm, m01,0);
	let (ngreater1, pval1, refdiff, nperm1) =
	    cluster_by_snp::test_snp_perm(&cl, &mmeth, &ncpgs, maxperm, m01,1);
	print!("{}\t{}\t{}", ngreater0, pval0, nperm0);
	println!("\t{}\t{}\t{}\t{}", ngreater1, pval1, nperm1, refdiff);
    }
}

fn filter_by_ncpgs<T:Copy>(v:&[T], ncpgs:&[usize], thr:usize) -> Vec<T>   {
    return (*v).iter().
    enumerate().
    filter(|x| ncpgs[(*x).0]> thr).
    map(|x| *(x.1)).collect();
}

fn main() {
    let ncpgs= [1,2,7,5,10,1,3,4,1,2];
    let mmeth = [0.4,0.5,0.6,0.2,0.1,0.3,0.4,0.5,0.6,0.7];
    let cl = [0,0,0,0,0,1,1,1,1,1];
    println!("{:?}", cl);
    let cl:Vec<usize> = filter_by_ncpgs(&cl, &ncpgs, 2); 
    println!("{:?}", cl);
    let mmeth:Vec<f32> =  filter_by_ncpgs(&mmeth, &ncpgs, 2); 
    println!("{:?}", mmeth);
    let ncpgs:Vec<usize> = filter_by_ncpgs(&ncpgs, &ncpgs, 2);
    println!("{:?}", ncpgs);
}

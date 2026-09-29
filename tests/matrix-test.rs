use asms::matrix;


#[test]
fn matrix_of_file_works() {
    let fname = "tests/resources/matrix-chr15-100863815-100863868.txt";
    let m = matrix::matrix_of_file(fname);
    assert!( m.idxgpos[1] == 2 );
    assert!( m.idxreads[1] == 0 );
    assert!( m.m[0*m.d+1] == 0.0 );
}

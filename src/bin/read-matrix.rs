extern crate asms;
use asms::matrix;

fn main() {
    let fnmatrix = "/home/ecmra/devel-ciacco/methylation/asms/paper/experiments/scan/chr15/matrices/matrix-chr15-100863815-100863868.txt";
    let mut _matrix;
    for _ in 0..10{
	_matrix = matrix::matrix_of_file(&fnmatrix);
    }
}

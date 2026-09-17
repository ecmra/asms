
use chrono;
use std::{env, io::Write, fs};
// see https://aeshirey.github.io/code/2021/02/13/adding-build-time-into-your-rust-package.html
fn main() {
    let outdir = env::var("OUT_DIR").unwrap();
    let outfile = format!("{}/timestamp.txt", outdir);
    let mut fh = fs::File::create(&outfile).unwrap();
    let current_date = chrono::Local::now().to_rfc2822();
    write!(fh, r#""{}""#, current_date).ok();
    //println!("cargo:rustc-link-lib=static=/lib/x86_64-linux-gnu/libssl.so.3")
}

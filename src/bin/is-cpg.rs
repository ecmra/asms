use std::io;

use std::io::BufRead;

use asms::utils;

fn main()  {
    let cpg =
	utils::read_cpg_list("/home/ecmra/devel-ciacco/methylation/solverd/data/human-cg.txt").
	unwrap();
    eprintln!("ready to read...");
    for line in io::stdin().lock().lines()  {
	match line {
	    Ok(line) => {
		let is = cpg.get(&line);
		match is {
		    Some(s) => println!("{}:{}", line, s),
		    _ => println!("{}:not found", line),
		}
	    },
	    Err(_e) => { eprintln!("we're done here!"); break}
	}
    }
    eprintln!("out of the loop");
}

//! find region at intermediate methylation

//! first scan: find all clusters of nadj consecutive
//! cpgs which are intermediately methylated.
//! second scan: merge such cluster when they are closer than
//! merge_distance

use crate::utils;

use rust_htslib::tbx;

use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::time;

#[derive(PartialEq)]
enum State {
    IN,
    OUT,
}

#[derive(PartialEq)]
enum Test {
    OK,
    NO,
}

use State::{IN, OUT};
use Test::{NO, OK};

pub fn scan(args: &[String]) {
    let usage = "usage: scan [OPTIONS] <methfn> <outdir>
<methfn>         compressed (.gz) tabix indexed methylbed file
<outdir>         directory which will contain the results
OPTIONS:
--lower-bound    <float in [0,1] default=0.4>
--upper-bound    <float in [0,1] default=0.6>
                 intermediate regions are defined as having methylation
                 between lower bound and upper bound.
--merge-distance <DISTANCE(bp)>
                 merge intermediate regions at less then this distance
                 (default=1000)
--nadj           <NCPGS>
                 minimum number of adjacent CpGs with intermediate methylation
                 to call an intermediate region.
                 (default=3)
--region         <CONTIG:NUM-NUM>
                 only look at specific region 
";
    if (0 == args.len()) || ("-h" == args[0]) || ("--help" == args[0]) {
        eprintln!("{}", usage);
        return;
    }
    let mut merge_distance: u32 = 1000;
    let mut nadj: usize = 3;
    let mut region: &str = "";
    let methfn: &String;
    let outdir: &String;
    let mut i: usize = 0;
    let mut lb:f32 = 0.4;
    let mut ub:f32 = 0.6;
    loop {
        if i >= args.len() {
            eprintln!("{}", usage);
            return;
        }
        if "--merge-distance" == args[i] {
            merge_distance = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
        if "--nadj" == args[i] {
            nadj = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
        if "--region" == args[i] {
            region = &args[i + 1];
            i = i + 2;
            continue;
        }
	if "--lower-bound" == args[i] {
            lb = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
	if "--upper-bound" == args[i] {
            ub = args[i + 1].parse().unwrap();
            i = i + 2;
            continue;
        }
        break;
    }
    if args.len() - i != 2 {
        eprintln!("{}", usage);
        return;
    }
    methfn = &args[i];
    outdir = &args[i + 1];
    eprintln!("nadj:{} ...", nadj);
    eprintln!("merge-distance:{}", merge_distance);
    let now = time::Instant::now();
    intermediate_regions(methfn, outdir, nadj, region, lb, ub);
    let imbedfn = format!("{}/intermediate-meth.bed", outdir);
    merge_nearby_im(&imbedfn, outdir, merge_distance);
    let elapsed = now.elapsed();
    eprintln!("scan ran in {:.3}sec(s)", elapsed.as_secs_f32());
}

/// find all clusters of <nadj> consecutive intermediately methylate CpGs.
fn intermediate_regions(methfn: &str, outdir: &str, nadj: usize, region: &str, lb:f32, ub:f32) -> () {
    eprintln!("trying to read from {} ...", methfn);
    eprintln!("writing results to {} ...", outdir);
    let mut tbxreader =
        tbx::Reader::from_path(&methfn).expect(&format!("Could not open {}", methfn));
    let outfn = format!("{}/intermediate-meth.bed", outdir);
    fs::create_dir_all(outdir).unwrap();
    let outfile = File::create(outfn).unwrap();
    let mut handle = BufWriter::new(outfile);
    let contigs = tbxreader.seqnames();
    let mut ms: utils::MethSlice;
    utils::write_out(&mut handle, &format!("#@NADJ:{}\n", nadj));
    if "" == region {
        for contig in contigs {
            eprintln!("Processing {}", contig);
            ms = utils::get_bedmethyl_slice(&mut tbxreader, &contig, 0, i64::MAX as u64);
            state_of_meth(&mut handle, ms, nadj, lb, ub);
        }
    } else {
        let (contig, start, end) = utils::parse_region(region);
        eprintln!("Processing {}:{}-{}", contig, start, end);
        ms = utils::get_bedmethyl_slice(&mut tbxreader, &contig, start as u64, end as u64);
        state_of_meth(&mut handle, ms, nadj, lb, ub);
    }
}

/// abstract the test
///
/// many possible tests specified by method
/// 0 default deterministic test
fn test_state(m: f32, lb: f32, ub: f32, method: u8) -> Test {
    match method {
        0 => match m {
            m if m >= lb && m <= ub => OK,
            _ => NO,
        },
        _ => panic!("not implemented"),
    }
}

/// generate first list of intermediate methylation regions
///
/// there must be at least nadj adjacent CpGs with intermediate
/// methylation to call a region.
fn state_of_meth<W: std::io::Write>(handle: &mut BufWriter<W>, ms: utils::MethSlice, nadj: usize, lb:f32, ub:f32) {
    let mut ncpg: usize = 0;
    let mut test: Test = NO;
    let mut state = OUT;
    let (mut irstart, mut irend) = (0u64, 0u64);
    let mut irlen: u64;
    for offset in 0..ms.meth.len() {
        test = test_state(ms.meth[offset], lb, ub, 0);
        state = match (&test, state) {
            (OK, IN) => {
                irend = ms.pos[offset];
                ncpg = ncpg + 1;
                IN
            }
            (OK, OUT) => {
                irstart = ms.pos[offset];
                irend = ms.pos[offset];
                ncpg = 1;
                IN
            }
            (NO, IN) => {
                irlen = irend - irstart + 1;
                if ncpg >= nadj && irlen >= 3 {
                    utils::write_out(
                        handle,
                        &format!(
                            "{}\t{}\t{}\t{}\t{}\n",
                            ms.contig,
                            irstart,
                            irend + 1,
                            irlen,
                            ncpg
                        ),
                    );
                }
                OUT
            }
            (NO, OUT) => OUT,
        };
    }
    if state == IN && test == OK {
        if ncpg >= nadj {
            irlen = irend - irstart + 1;
            utils::write_out(
                handle,
                &format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    ms.contig,
                    irstart,
                    irend + 1,
                    irlen,
                    ncpg
                ),
            );
        }
    }
}

fn merge_nearby_im(imbed: &str, outdir: &str, thr: u32) {
    let outfn = format!("{}/merged-im.bed", outdir);
    let outfile = File::create(outfn).unwrap();
    let mut handle = BufWriter::new(outfile);
    utils::write_out(&mut handle, &format!("#@MERGE-DISTANCE:{}\n", thr));
    let mut bufreader = BufReader::new(File::open(imbed).unwrap());
    eprintln!("reading from {} ...", imbed);
    let mut buffer = String::new();
    let mut len: usize;
    let mut res;
    loop {
        res = bufreader.read_line(&mut buffer);
        len = match res {
            Ok(l) => l,
            _ => panic!("problem in reading bed file"),
        };
        if 0 == len {
            panic!("empty bed file");
        }
        if buffer.starts_with("#") {
            buffer.clear();
            continue;
        } else {
            break;
        }
    }
    let oldline = buffer.clone();
    let mut fields: Vec<&str> = vec![""; 5];
    let mut oldfields = oldline
        .trim_end()
        .split_whitespace()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    let mut line: String;
    let mut state = OUT;
    buffer.clear();
    let mut left: u32 = 0;
    let mut right: u32;
    let mut cpgcount: u32 = 0;

    loop {
        res = bufreader.read_line(&mut buffer);
        len = match res {
            Ok(l) => l,
            _ => panic!("error in reading bed file"),
        };
        if 0 == len {
            break;
        }
        line = buffer.clone();
        buffer.clear();
        fields = line.trim_end().split_whitespace().collect::<Vec<_>>();
        let dist = if fields[0] == oldfields[0] {
            fields[1].parse::<u32>().unwrap() - oldfields[2].parse::<u32>().unwrap()
        } else {
            thr + 1
        };
        state = match (dist, state) {
            (d, s) if d <= thr && s == OUT => {
                left = oldfields[1].parse::<u32>().unwrap();
                cpgcount = oldfields[4].parse::<u32>().unwrap() + fields[4].parse::<u32>().unwrap();
                IN
            }
            (d, s) if d <= thr && s == IN => {
                cpgcount = cpgcount + fields[4].parse::<u32>().unwrap();
                IN
            }
            (d, s) if d > thr && s == OUT => {
                left = oldfields[1].parse::<u32>().unwrap();
                right = oldfields[2].parse::<u32>().unwrap();
                cpgcount = oldfields[4].parse::<u32>().unwrap();
                utils::write_out(
                    &mut handle,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\n",
                        oldfields[0],
                        left,
                        right,
                        right - left,
                        cpgcount
                    ),
                );
                OUT
            }
            (d, s) if d > thr && s == IN => {
                right = oldfields[2].parse::<u32>().unwrap();
                utils::write_out(
                    &mut handle,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\n",
                        oldfields[0],
                        left,
                        right,
                        right - left,
                        cpgcount
                    ),
                );
                OUT
            }
            (_, _) => panic!("uh?"),
        };
        for j in 0..fields.len() {
            oldfields[j] = String::from(fields[j]);
        }
    }
    // on exit
    eprintln!("exit");
    match state {
        IN => {
            right = fields[2].parse::<u32>().unwrap();
            utils::write_out(
                &mut handle,
                &format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    oldfields[0],
                    left,
                    right,
                    right - left,
                    cpgcount
                ),
            );
        }
        OUT => {
            left = fields[1].parse::<u32>().unwrap();
            right = fields[2].parse::<u32>().unwrap();
            cpgcount = fields[4].parse::<u32>().unwrap();
            utils::write_out(
                &mut handle,
                &format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    oldfields[0],
                    left,
                    right,
                    right - left,
                    cpgcount
                ),
            );
        }
    };
}

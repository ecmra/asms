use rust_htslib::htslib;
use rust_htslib::bam::record::Seq;

/// base modification iter
/// can give positions which correspond to a G on the read
/// when the mod is detected on the other strand.
pub fn parse_mods(record: &rust_htslib::bam::Record,
	      calls: &mut Vec<(i32, htslib::hts_base_mod)>, seq:Seq<'_>) {
    let rname = std::str::from_utf8(record.qname()).unwrap();
    /*     mods has all the methods of an iterator,
           with item = Result<(i32 (pos), hts_base_mod)>
           plus extra methods to extract info
           about how many and what type of
           modifications are recorded fot the record
           I am looking at
     */
    let mods = match record.basemods_iter() {
        Ok(mods) => mods,
        Err(_) => {
            println!("Failed to parse MM/ML aux tags in {}\n", rname);
            return;
        }
    };
    let mcodes: &[i32] = mods.recorded();
    let mut implicit: bool = true;
    let mut canbase: u8 = b'C';
    /* there can be many modifications on the same position.
    we look for 'm' and keep track of whether
    it is stored as implicit and what the canonical base is */
    for i in 0..mcodes.len() {
        let modtype = mods.query_type(mcodes[i]).unwrap();
        canbase = modtype.canonical;
        if (canbase == b'C' || canbase == b'G')
            && mcodes[i] as u8  == b'm'
            && modtype.implicit == 0
        {
            implicit = false;
            break;
        };
    }
    let mut explpos: Vec<i32> = vec![];
    for res in mods {
        if let Ok((position, m)) = res {
            if m.modified_base as u8  != b'm' {
                continue;
            }
            calls.push((position, m));
            explpos.push(position);
        }
    }
    /* add implicit calls: scan the sequence in the 5' to 3' direction.
       but add positions in the coordinate of the SEQ
       as it appears in the BAM file.
       Add the C/G corresponding to implicit
       absence of methylation. Used when methylation
       is encoded implicitly, eg C+m. or G-m. */
    if implicit {
        let isrev = record.flags() & 0x10 == 0x10;
	//let seq = seq.as_bytes();
        let target = match (canbase, isrev) {
            (b'C', false) => b'C',
            (b'C', true)  => b'G',
            (b'G', false) => b'G',
            (b'G', true)  => b'C',
            _ => panic!("unrecognized canonical base {}{:?}",
				canbase, isrev),
        };
        for pos in 0..record.seq_len() {
            if unsafe{ seq.decoded_base_unchecked(pos) }  == target  {
                if !(explpos.contains(&(pos as i32))) {
                    let mb = htslib::hts_base_mod {
                        modified_base: b'm' as i32,
                        canonical_base: canbase as i32,
                        strand: match canbase  {
                            b'C' => 0,
                            b'G' => 1,
                            _ => panic!("fill_implicit::panic!"),
                        },
                        qual: 0,
                    };
                    calls.push((pos as i32, mb));
                }
            }
        }
    }
}

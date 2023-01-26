fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::io::*;
    env_logger::init();
    let mut diff = vec![];
    let records: Vec<_> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|l| bio_utils::paf::PAF::new(&l))
        .collect();
    for record in records.iter() {
        let mut qpos = 0;
        let mut rpos = record.tstart;
        let qname = &record.qname;
        let rname = &record.tname;
        let cigar = bio_utils::sam::parse_cigar_string(record.get_tag("cg").unwrap().1);
        for op in cigar {
            use bio_utils::sam::Op;
            match op {
                Op::Match(l) | Op::Align(l) => {
                    qpos += l;
                    rpos += l;
                }
                Op::Insertion(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Ins", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Ins", l));
                    }
                    qpos += l;
                }
                Op::Deletion(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Del", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Del", l));
                    }
                    rpos += l;
                }
                Op::Mismatch(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Mism", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Mism", l));
                    }
                    rpos += l;
                    qpos += l;
                }
                _ => {}
            }
        }
    }
    println!("Query\tQpos\tRefr\tRpos\tType\tSize");
    for (qname, qpos, rname, rpos, t, size) in diff.iter() {
        println!("{qname}\t{qpos}\t{rname}\t{rpos}\t{t}\t{size}");
    }
    Ok(())
}

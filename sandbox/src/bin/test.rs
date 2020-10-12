fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let sam: Vec<_> = bio_utils::sam::load_sam_file(&args[1]).unwrap();
    // let _reference: HashMap<_, _> = bio_utils::fasta::parse_into_vec(&args[2])
    //     .unwrap()
    //     .into_iter()
    //     .map(|r| (r.id().to_string(), r))
    //     .collect();
    println!("RefName\tMat\tMism\tIns\tDel");
    for record in sam {
        let (mut ins, mut del, mut mat, mut mism) = (0, 0, 0, 0);
        use bio_utils::sam::Op;
        for op in record.cigar() {
            match op {
                Op::Insertion(l) => ins += l,
                Op::Deletion(l) => del += l,
                Op::Match(l) => mat += l,
                Op::Align(l) => mat += l,
                Op::Mismatch(l) => mism += l,
                _ => {}
            }
        }
        let (mut max_ins, mut max_del) = (0, 0);
        let (mut r_pos, mut q_pos) = (record.pos(), 0);
        for op in record.cigar() {
            match op {
                Op::Insertion(l) => {
                    if l >= 80 {
                        eprintln!("{}\t{}\t{}\tIns", q_pos, r_pos, l);
                    }
                    max_ins = max_ins.max(l);
                    q_pos += l;
                }
                Op::Deletion(l) => {
                    if l >= 80 {
                        eprintln!("{}\t{}\t{}\tDel", q_pos, r_pos, l);
                    }
                    max_del = max_del.max(l);
                    r_pos += l;
                }
                Op::Match(l) | Op::Align(l) | Op::Mismatch(l) => {
                    r_pos += l;
                    q_pos += l;
                }
                Op::SoftClip(l) => q_pos += l,
                _ => {}
            }
        }
        // println!("{}\t{}", max_ins, max_del);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            record.r_name(),
            record.q_name(),
            mat,
            mism,
            ins,
            del
        );
        // let seq = record.seq().as_bytes();
        // if seq.len() > 1 {
        //     let cigar = record.cigar();
        //     let refr = reference[record.r_name()].seq();
        //     let (q, o, r) = bio_utils::sam::recover_alignment(&cigar, seq, &refr, record.pos());
        //     for ((q, o), r) in q.chunks(150).zip(o.chunks(150)).zip(r.chunks(150)) {
        //         println!("{}", String::from_utf8_lossy(q));
        //         println!("{}", String::from_utf8_lossy(o));
        //         println!("{}\n", String::from_utf8_lossy(r));
        //     }
        // }
    }
}

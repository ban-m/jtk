fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();

    use definitions::*;
    use std::io::*;
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let target = 274;
    let nodes: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|r| r.nodes.iter())
        .filter(|n| n.unit == target)
        .collect();
    let mut hmm = haplotyper::model_tune::get_model(&ds).unwrap();
    let mut cons = ds
        .selected_chunks
        .iter()
        .find(|u| u.id == target)
        .unwrap()
        .seq()
        .to_vec();
    let mut seqs: Vec<_> = nodes.iter().map(|n| n.seq().to_vec()).collect();
    let mut ops: Vec<_> = nodes
        .iter()
        .map(|n| haplotyper::misc::ops_to_kiley(&n.cigar))
        .collect();
    let band_width = 100;
    let total = seqs.iter().map(|x| x.len()).sum::<usize>();
    eprintln!("{total}");
    for _t in 0..2 {
        for (seq, ops) in seqs.iter_mut().zip(ops.iter_mut()) {
            haplotyper::misc::fix_long_homopolymers(seq, &cons, ops, 5);
        }
        hmm.fit_naive_with(&cons, &seqs, &ops, band_width / 2);
        cons = hmm.polish_until_converge_with(&cons, &seqs, &mut ops, band_width / 2);
    }
    let total = seqs.iter().map(|x| x.len()).sum::<usize>();
    eprintln!("{total}");
    // for (seq, ops) in seqs.iter().zip(ops.iter()) {
    // let (xr, ar, yr) = kiley::recover(&cons, &seq, &ops);
    // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //     eprintln!("{}", std::str::from_utf8(xr).unwrap());
    //     eprintln!("{}", std::str::from_utf8(ar).unwrap());
    //     eprintln!("{}\n", std::str::from_utf8(yr).unwrap());
    // }
    // eprintln!("==================")
    //    }
    for (i, seq) in seqs.iter().enumerate() {
        println!("READ\t{i}\t{}", std::str::from_utf8(seq).unwrap());
    }
    println!("UNIT\t{target}\t{}", std::str::from_utf8(&cons).unwrap());
    // use std::collections::HashSet;
    // let selection: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    // use haplotyper::local_clustering::*;
    // local_clustering_selected(&mut ds, &selection);
    // println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // let mut diff = vec![];
    // let records: Vec<_> = std::fs::File::open(&args[1])
    //     .map(BufReader::new)?
    //     .lines()
    //     .filter_map(|x| x.ok())
    //     .filter_map(|l| bio_utils::paf::PAF::new(&l))
    //     .collect();
    // for record in records.iter() {
    //     let mut qpos = record.qstart;
    //     let mut rpos = record.tstart;
    //     let qname = &record.qname;
    //     let rname = &record.tname;
    //     let cigar = bio_utils::sam::parse_cigar_string(&record.get_tag("cg").unwrap().1);
    //     for op in cigar {
    //         use bio_utils::sam::Op;
    //         match op {
    //             Op::Match(l) | Op::Align(l) => {
    //                 qpos += l;
    //                 rpos += l;
    //             }
    //             Op::Insertion(l) => {
    //                 if record.relstrand {
    //                     diff.push((qname, qpos, rname, rpos, "Ins", l));
    //                 } else {
    //                     let qpos = record.qlen - qpos;
    //                     diff.push((qname, qpos, rname, rpos, "Ins", l));
    //                 }
    //                 qpos += l;
    //             }
    //             Op::Deletion(l) => {
    //                 if record.relstrand {
    //                     diff.push((qname, qpos, rname, rpos, "Del", l));
    //                 } else {
    //                     let qpos = record.qlen - qpos;
    //                     diff.push((qname, qpos, rname, rpos, "Del", l));
    //                 }
    //                 rpos += l;
    //             }
    //             Op::Mismatch(l) => {
    //                 if record.relstrand {
    //                     diff.push((qname, qpos, rname, rpos, "Mism", l));
    //                 } else {
    //                     let qpos = record.qlen - qpos;
    //                     diff.push((qname, qpos, rname, rpos, "Mism", l));
    //                 }
    //                 rpos += l;
    //                 qpos += l;
    //             }
    //             _ => {}
    //         }
    //     }
    // }
    // println!("Query\tQpos\tRefr\tRpos\tType\tSize");
    // for (qname, qpos, rname, rpos, t, size) in diff.iter() {
    //     println!("{qname}\t{qpos}\t{rname}\t{rpos}\t{t}\t{size}");
    // }
    Ok(())
}

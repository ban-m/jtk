fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    use std::fs::File;
    let mut wtr_a = File::create("hapA.fa").map(bio_utils::fasta::Writer::new)?;
    let mut wtr_c = File::create("hapC.fa").map(bio_utils::fasta::Writer::new)?;
    for read in reads {
        if read.id().contains("hapA") {
            wtr_a.write_record(&read)?;
        } else {
            wtr_c.write_record(&read)?;
        }
    }
    // let mut reference = bio_utils::fasta::parse_into_vec(&args[1])?;
    // let reads = bio_utils::fasta::parse_into_vec(&args[2])?;
    // let reads: Vec<_> = reads
    //     .into_iter()
    //     .enumerate()
    //     .map(|(idx, read)| definitions::RawRead {
    //         name: read.id().to_string(),
    //         desc: read.desc().map(|i| i.to_string()).unwrap_or(String::new()),
    //         id: idx as u64,
    //         seq: String::from_utf8_lossy(read.seq()).to_string(),
    //     })
    //     .collect();
    // let mut cdir = std::env::current_dir().unwrap();
    // cdir.push("maptemp");
    // std::fs::create_dir_all(&cdir)?;
    // let alignments: Vec<_> = haplotyper::assemble::invoke_last(&args[1], &args[2], &cdir, 24)?;
    // let flag = true;
    // if flag {
    //     let segment = {
    //         let refr = reference.pop().unwrap();
    //         let seq = Some(String::from_utf8_lossy(refr.seq()).to_string());
    //         gfa::Segment::from(refr.id().to_string(), refr.seq().len(), seq)
    //     };
    //     let reads: Vec<_> = reads.iter().collect();
    //     let c = haplotyper::assemble::AssembleConfig::new(24, 100);
    //     let result = haplotyper::assemble::polish_by_chunking(&alignments, &segment, &reads, &c);
    //     let result = bio_utils::fasta::Record::with_data("tig", &None, &result);
    //     println!("{}", result);
    // } else {
    //     const LOCATION: usize = 809_699;
    //     use std::collections::HashMap;
    //     let mut bucket: HashMap<_, Vec<_>> = HashMap::new();
    //     for aln in alignments.iter() {
    //         let key = aln.seq2_name().to_string();
    //         bucket.entry(key).or_default().push(aln);
    //     }
    //     let references: HashMap<_, _> = reference
    //         .iter()
    //         .map(|r| (r.id().to_string(), r.seq()))
    //         .collect();
    //     use haplotyper::encode::join_alignments;
    //     for read in reads.iter() {
    //         if let Some(alns) = bucket.get(&read.name) {
    //             if let Some(refr) = references.get(alns[0].seq1_name()) {
    //                 for &dir in &[false, true] {
    //                     let alns: Vec<_> = alns
    //                         .iter()
    //                         .copied()
    //                         .filter(|e| e.seq2_direction().is_forward() == dir)
    //                         .collect();
    //                     let seq = if dir {
    //                         read.seq().to_vec()
    //                     } else {
    //                         bio_utils::revcmp(read.seq())
    //                     };
    //                     if !alns.is_empty() {
    //                         use definitions::Op;
    //                         let (qpos, rpos, ops) = join_alignments(&alns, refr, &seq);
    //                         let matchlen = ops
    //                             .iter()
    //                             .map(|op| match *op {
    //                                 Op::Match(l) | Op::Del(l) => l,
    //                                 Op::Ins(_) => 0,
    //                             })
    //                             .sum::<usize>();
    //                         if rpos < LOCATION && LOCATION < rpos + matchlen {
    //                             for aln in alns.iter() {
    //                                 println!("{}", aln);
    //                             }
    //                             println!("{}-{}", rpos, rpos + matchlen);
    //                             // let (r, o, q) =
    //                             //     haplotyper::encode::recover(&seq[qpos..], &refr[rpos..], &ops);
    //                             let chunks =
    //                                 haplotyper::assemble::split_reads(&seq, qpos, rpos, ops, 100);
    //                             for (idx, c) in chunks.iter() {
    //                                 println!("{}\t{}", idx, c.len());
    //                             }
    //                             // for ((r, o), q) in r.chunks(150).zip(o.chunks(150)).zip(q.chunks(150)) {
    //                             //     println!("{}", String::from_utf8_lossy(r));
    //                             //     println!("{}", String::from_utf8_lossy(o));
    //                             //     println!("{}", String::from_utf8_lossy(q));
    //                             //     println!();
    //                             // }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    Ok(())
}

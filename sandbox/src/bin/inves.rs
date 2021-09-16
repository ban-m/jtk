use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    let target = 376;
    let ref_unit = ds.selected_chunks.iter().find(|u| u.id == target).unwrap();
    use std::collections::HashMap;
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|r| r.nodes.iter())
        .filter(|n| n.unit == target)
        .map(|n| n.seq())
        .collect();
    let band_width = 100;
    let cons_template = kiley::consensus(&reads, 1000, 10, band_width);
    println!(
        "DUP:{}",
        String::from_utf8_lossy(&cons_template[1183..1187])
    );
    eprintln!(
        "SEQ\t>prv\nSEQ\t{}",
        String::from_utf8_lossy(&cons_template)
    );
    let dup_template: Vec<_> = cons_template[..1187]
        .iter()
        .chain(cons_template[1183..].iter())
        .copied()
        .collect();
    eprintln!("SEQ\t>aft\nSEQ\t{}", String::from_utf8_lossy(&dup_template));
    for read in reads {
        println!();
        let (score, ops) =
            kiley::bialignment::global_banded(&cons_template, read, 2, -5, -6, -1, 100);
        let (xr, ar, yr) = kiley::bialignment::recover(&cons_template, read, &ops);
        let (mut cbuf, mut abuf, mut rbuf) = (String::new(), String::new(), String::new());
        let mut cpos = 0;
        for ((&xr, &ar), &yr) in xr.iter().zip(ar.iter()).zip(yr.iter()) {
            cpos += (xr != b' ') as usize;
            if (1170..1270).contains(&cpos) {
                cbuf.push(xr as char);
                abuf.push(ar as char);
                rbuf.push(yr as char);
            }
        }
        println!("{}\n{}\n{}\n{}", score, cbuf, abuf, rbuf);
        let (score, ops) =
            kiley::bialignment::global_banded(&dup_template, read, 2, -5, -6, -1, 100);
        let (xr, ar, yr) = kiley::bialignment::recover(&dup_template, read, &ops);
        let (mut cbuf, mut abuf, mut rbuf) = (String::new(), String::new(), String::new());
        let mut cpos = 0;
        for ((&xr, &ar), &yr) in xr.iter().zip(ar.iter()).zip(yr.iter()) {
            cpos += (xr != b' ') as usize;
            if (1170..1270).contains(&cpos) {
                cbuf.push(xr as char);
                abuf.push(ar as char);
                rbuf.push(yr as char);
            }
        }
        println!("{}\n{}\n{}\n{}\n", score, cbuf, abuf, rbuf);
    }
    // for idx in [25, 35, 1183, 1184, 1185, 1186] {
    //     for rep_size in 1..4 {
    //         let mod_template: Vec<_> = cons_template[..idx + rep_size]
    //             .iter()
    //             .chain(cons_template[idx..].iter())
    //             .copied()
    //             .collect();
    //         for (i, read) in reads.iter().enumerate() {
    //             let (prev, _) =
    //                 kiley::bialignment::global_banded(&cons_template, read, 2, -5, -6, -1, 100);
    //             let (after, _) =
    //                 kiley::bialignment::global_banded(&mod_template, read, 2, -5, -6, -1, 100);
    //             println!("{}\t{}\t{}\t{}", i, idx, rep_size, after - prev);
    //         }
    //     }
    // }
    // use kiley::gphmm::*;
    // let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    // let template = kiley::padseq::PadSeq::new(cons_template.as_slice());
    // let reads: Vec<_> = reads
    //     .iter()
    //     .map(|&r| kiley::padseq::PadSeq::new(r))
    //     .collect();
    // let (hmm, _) = hmm.fit_banded_inner(&template, &reads, band_width as usize);
    // const REP_SIZE: usize = 5;
    // let profiles: Vec<_> = reads
    //     .iter()
    //     .map(|read| {
    //         let prof =
    //             banded::ProfileBanded::new(&hmm, &template, read, band_width as isize).unwrap();
    //         prof.to_copy_table(REP_SIZE)
    //     })
    //     .collect();
    // for pos in 0..REP_SIZE * (template.len() - REP_SIZE) {
    //     let rep_size = pos % REP_SIZE + 1;
    //     let idx = pos / REP_SIZE;
    //     if ![25, 35, 1183, 1184, 1185, 1186].contains(&idx) {
    //         continue;
    //     }
    //     println!("{},{},{}", pos, idx, rep_size);
    //     // So, modify in such way.
    //     let mod_template: Vec<_> = cons_template[..idx + rep_size]
    //         .iter()
    //         .chain(cons_template[idx..].iter())
    //         .copied()
    //         .collect();
    //     let mod_template = kiley::padseq::PadSeq::new(mod_template.as_slice());
    //     for (read, prof) in reads.iter().zip(profiles.iter()) {
    //         assert_eq!(prof.len(), REP_SIZE * (template.len() - REP_SIZE));
    //         let mod_lk = hmm.likelihood_inner(&mod_template, &read);
    //         let diff = (mod_lk - prof[pos]).abs();
    //         assert!(diff < 0.001, "{},{}", mod_lk, prof[pos]);
    //     }
    // }
    // let (hap1, hap2): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .filter(|read| read.nodes.iter().any(|n| n.unit == target))
    //     .partition(|read| id2desc[&read.id].contains("000252v2"));
    // let reads = hap1
    //     .iter()
    //     .map(|r| (0, r))
    //     .chain(hap2.iter().map(|r| (1, r)));
    // for (hap, read) in reads {
    //     for node in read.nodes.iter().filter(|n| n.unit == target) {
    //         println!("======={}========", hap);
    //         let (xr, ar, yr) = node.recover(ref_unit);
    //         for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //             println!("{}", String::from_utf8_lossy(xr));
    //             println!("{}", String::from_utf8_lossy(ar));
    //             println!("{}\n", String::from_utf8_lossy(yr));
    //         }
    // println!("===========POLISH============");
    // let (_, ops) =
    //     kiley::bialignment::global_banded(node.seq(), ref_unit.seq(), 2, -6, -5, -1, 100);
    // let (xr, ar, yr) = kiley::bialignment::recover(node.seq(), ref_unit.seq(), &ops);
    // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //     println!("{}", String::from_utf8_lossy(xr));
    //     println!("{}", String::from_utf8_lossy(ar));
    //     println!("{}\n", String::from_utf8_lossy(yr));
    // }
    //     }
    // }
    // use haplotyper::assemble::*;
    // let config = AssembleConfig::new(3, 1000, false, false, 4);
    // ds.squish_small_contig(&config, 30);
    // println!("{}", ds.assemble(&config));
    // let counts: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter_map(|r| {
    //         let nodes: Vec<_> = r
    //             .nodes
    //             .iter()
    //             .filter_map(|n| [905, 924].contains(&n.unit).then(|| (n.unit, n.cluster)))
    //             .collect();
    //         (nodes.len() == 2).then(|| nodes)
    //     })
    //     .collect();
    // for nodes in counts {
    //     println!("{:?}", nodes);
    // }

    // let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    // use haplotyper::assemble::copy_number::CoverageCalibrator;
    // let calib = CoverageCalibrator::new(&lens);
    // println!("gaplen\tspan");
    // for len in (0..).map(|i| 1_000 + 1_000 * i).take_while(|&x| x < 50_000) {
    //     println!("{}\t{}", len, calib.prob_spanning(len));
    // }
    // use kiley::sam::Sam;
    // let records = std::fs::File::open(&args[1])
    //     .map(BufReader::new)
    //     .map(Sam::from_reader)?;
    // let record = kiley::fasta::read_fasta(&Some(&args[2]))?;
    // let (mut xpos, mut ypos) = (0, 0);
    // let (xs, ys) = (&record[0].1, &record[1].1);
    // let ofs = 5;
    // for op in records.records[0].cigar().iter() {
    //     match op {
    //         kiley::sam::Op::Match(l) | kiley::sam::Op::Align(l) => {
    //             xpos += l;
    //             ypos += l;
    //         }
    //         kiley::sam::Op::Insertion(l) => {
    //             ypos += l;
    //             let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
    //             let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
    //             println!("{},{},{},I,{},{}", xpos, ypos, l, xseq, yseq);
    //         }
    //         kiley::sam::Op::Deletion(l) => {
    //             xpos += l;
    //             let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
    //             let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
    //             println!("{},{},{},D,{},{}", xpos, ypos, l, xseq, yseq);
    //         }
    //         kiley::sam::Op::Mismatch(l) => {
    //             xpos += l;
    //             ypos += l;
    //             let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
    //             let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
    //             println!("{},{},{},M,{},{}", xpos, ypos, l, xseq, yseq);
    //         }
    //         kiley::sam::Op::SoftClip(l) | kiley::sam::Op::HardClip(l) => {
    //             ypos += l;
    //         }
    //         _ => {}
    //     }
    // }
    Ok(())
}

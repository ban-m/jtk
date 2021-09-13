use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    use haplotyper::assemble::*;
    let config = AssembleConfig::new(3, 1000, false, true, 4);
    // ds.squish_small_contig(&config, 30);
    println!("{}", ds.assemble(&config));
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

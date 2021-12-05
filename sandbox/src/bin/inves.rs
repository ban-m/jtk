use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // Assembly
    use haplotyper::assemble::*;
    let config = AssembleConfig::new(3, 1000, false, true, 4);
    println!("{}", ds.assemble(&config));
    // use std::io::BufRead;
    // let input: Vec<_> = std::fs::File::open(&args[1])
    //     .map(BufReader::new)?
    //     .lines()
    //     .filter_map(|l| l.ok())
    //     .collect();
    // let mut paf: HashMap<_, Vec<_>> = HashMap::new();
    // for line in input {
    //     let id: u64 = line.split('\t').next().unwrap().parse().unwrap();
    //     paf.entry(id).or_default().push(line);
    // }
    // paf.retain(|_, val| val.len() == 2);
    // let record_to_mntag = |record: &String| -> Option<i64> {
    //     record
    //         .split('\t')
    //         .find(|tag| tag.starts_with("NM"))?
    //         .rsplit_once(':')?
    //         .1
    //         .parse()
    //         .ok()
    // };
    // paf.retain(|_, val| {
    //     let nm: Vec<_> = val.iter().filter_map(record_to_mntag).collect();
    //     !nm.iter().all(|&x| x == nm[0])
    // });
    // let record_to_mism = |record: &String| -> Option<usize> {
    //     let cigar = record
    //         .split('\t')
    //         .find(|tag| tag.starts_with("cg"))?
    //         .rsplit_once(':')?
    //         .1;
    //     let mismatch: usize = bio_utils::sam::parse_cigar_string(&cigar)
    //         .iter()
    //         .filter_map(|op| match op {
    //             &bio_utils::sam::Op::Mismatch(x) => Some(x),
    //             _ => None,
    //         })
    //         .sum();
    //     Some(mismatch)
    // };
    // paf.retain(|_, val| {
    //     let mism: Vec<_> = val.iter().filter_map(record_to_mism).collect();
    //     !mism.iter().all(|&x| x == mism[0])
    // });
    // let mut paf: Vec<_> = paf.into_iter().collect();
    // paf.sort_by_key(|x| x.0);
    // for (_, bucket) in paf.iter() {
    //     for rec in bucket {
    //         println!("{}", rec);
    //     }
    // }
    // use definitions::*;
    // let ds: DataSet = std::fs::File::open(&args[1])
    //     .map(BufReader::new)
    //     .map(|x| serde_json::de::from_reader(x).unwrap())?;
    // let target: u64 = args[2].parse().unwrap();
    // let ref_unit = ds.selected_chunks.iter().find(|u| u.id == target).unwrap();
    // //let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // let (hap1, hap2): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .filter(|read| read.nodes.iter().any(|n| n.unit == target))
    //     .partition(|read| id2name[&read.id].contains("ler"));
    // // .partition(|read| id2desc[&read.id].contains("000252v2"));
    // let reads = hap1
    //     .iter()
    //     .map(|r| (0, r))
    //     .chain(hap2.iter().map(|r| (1, r)));
    // for (hap, read) in reads {
    //     for node in read.nodes.iter().filter(|n| n.unit == target) {
    //         let indel_iter = node.cigar.iter().map(|&op| match op {
    //             Op::Del(l) | Op::Ins(l) => l as i32,
    //             Op::Match(l) => -(l as i32),
    //         });
    //         use haplotyper::encode::max_region;
    //         let max_indel = max_region(indel_iter);
    //         let (mut npos, mut rpos) = (0, 0);
    //         let (nodeseq, refseq) = (node.seq(), ref_unit.seq());
    //         let mut mism = 0;
    //         for op in node.cigar.iter() {
    //             match *op {
    //                 Op::Match(l) => {
    //                     mism += nodeseq
    //                         .iter()
    //                         .skip(npos)
    //                         .take(l)
    //                         .zip(refseq.iter().skip(rpos).take(l))
    //                         .filter(|(x, y)| x != y)
    //                         .count();
    //                     rpos += l;
    //                     npos += l;
    //                 }
    //                 Op::Del(l) => rpos += l,
    //                 Op::Ins(l) => npos += l,
    //             }
    //         }
    //         //let desc = id2desc[&read.id];
    //         let desc = id2name[&read.id];
    //         let id = read.id;
    //         let pos = node.position_from_start;
    //         println!(
    //             "{},{},{},{},{},{},{}",
    //             hap, node.cluster, max_indel, mism, id, pos, desc
    //         );
    //         let (xr, ar, yr) = node.recover(ref_unit);
    //         for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //             println!("{}", String::from_utf8_lossy(xr));
    //             println!("{}", String::from_utf8_lossy(ar));
    //             println!("{}\n", String::from_utf8_lossy(yr));
    //         }
    // println!("===========POLISH============");
    // let (_, ops) =
    //     kiley::bialignment::global_banded(node.seq(), ref_unit.seq(), 2, -7, -4, -1, 100);
    // let indel_iter = ops.iter().map(|op| match op {
    //     kiley::bialignment::Op::Del => 1,
    //     kiley::bialignment::Op::Ins => 1,
    //     kiley::bialignment::Op::Mat => -1,
    // });
    // use haplotyper::encode::max_region;
    // let max_indel = max_region(indel_iter);
    // println!("======={},{},{}========", hap, node.cluster, max_indel);
    // let (xr, ar, yr) = kiley::bialignment::recover(node.seq(), ref_unit.seq(), &ops);
    // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //     println!("{}", String::from_utf8_lossy(xr));
    //     println!("{}", String::from_utf8_lossy(ar));
    //     println!("{}\n", String::from_utf8_lossy(yr));
    // }
    // let prev = node
    //     .cigar
    //     .iter()
    //     .map(|x| match x {
    //         Op::Match(_) => 0,
    //         Op::Del(l) | Op::Ins(l) => *l,
    //     })
    //     .max()
    //     .unwrap();
    // let (_, ops) =
    //     kiley::bialignment::global_banded(ref_unit.seq(), node.seq(), 2, -5, -6, -1, 200);
    // let aln_dist = {
    //     let (mut upos, mut spos) = (0, 0);
    //     ops.iter()
    //         .map(|op| match op {
    //             kiley::bialignment::Op::Del => {
    //                 upos += 1;
    //                 1
    //             }
    //             kiley::bialignment::Op::Ins => {
    //                 spos += 1;
    //                 1
    //             }
    //             kiley::bialignment::Op::Mat => {
    //                 spos += 1;
    //                 upos += 1;
    //                 (ref_unit.seq()[upos - 1] != node.seq()[spos - 1]) as u32
    //             }
    //         })
    //         .sum::<u32>()
    // };
    // let cigar = haplotyper::encode::compress_kiley_ops(&ops);
    // let after = cigar
    //     .iter()
    //     .map(|x| match x {
    //         Op::Match(_) => 0,
    //         Op::Del(l) | Op::Ins(l) => *l,
    //     })
    //     .max()
    //     .unwrap();
    // println!(
    //     "{}\t{}\t{}\t{}\t{}",
    //     hap, node.cluster, prev, after, aln_dist
    // );
    //     }
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

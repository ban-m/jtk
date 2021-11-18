use definitions::*;
use haplotyper::Encode;
// use haplotyper::DetermineUnit;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    use haplotyper::encode::deletion_fill;
    ds = deletion_fill::correct_unit_deletion(ds, 0.35);
    use haplotyper::re_clustering::*;
    let config = ReClusteringConfig::new(24, 7, 5);
    ds = ds.re_clustering(&config);
    use std::io::{BufWriter, Write};
    if let Ok(mut wtr) = std::fs::File::create("bf.json").map(BufWriter::new) {
        writeln!(
            &mut wtr,
            "{}",
            serde_json::ser::to_string_pretty(&ds).unwrap()
        )?;
    }
    use haplotyper::dense_encoding::*;
    let config = DenseEncodingConfig::new(5, 5);
    ds = ds.dense_encoding_hapcut(&config);
    if let Ok(mut wtr) = std::fs::File::create("af.json").map(BufWriter::new) {
        writeln!(
            &mut wtr,
            "{}",
            serde_json::ser::to_string_pretty(&ds).unwrap()
        )?;
    }
    ds = deletion_fill::correct_unit_deletion(ds, 0.35);
    println!("{}", serde_json::ser::to_string_pretty(&ds).unwrap());
    // // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // for read in ds.encoded_reads.iter_mut() {
    //     for node in read.nodes.iter_mut() {
    //         node.cluster = 0;
    //     }
    // }
    // use assemble::Assemble;
    // use haplotyper::assemble;
    // let config = assemble::AssembleConfig::new(10, 1000, false, false, 5);
    // let gfa = ds.assemble(&config);
    // println!("{}", gfa);
    // let units: Vec<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    // for read in ds.encoded_reads.iter() {
    //     if units
    //         .iter()
    //         .all(|&n| read.nodes.iter().any(|node| node.unit == n))
    //     {
    //         let answer = id2desc[&read.id].contains("251v2") as usize;
    //         let read: Vec<_> = read
    //             .nodes
    //             .iter()
    //             .filter(|n| units.contains(&n.unit))
    //             .map(|n| format!("{}-{}", n.unit, n.cluster))
    //             .collect();
    //         println!("{}\t{}", answer, read.join("\t"));
    //     }
    // }
    // Length of the reads given.
    // Match names.
    // Length of homopolymers.
    // let mut homo_polymers: Vec<Vec<u32>> = vec![vec![]; 4];
    // for seq in ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| r.seq())
    //     .filter(|x| !x.is_empty())
    // {
    //     let mut seq = seq.iter().map(|x| x.to_ascii_uppercase()).peekable();
    //     while let Some(current_base) = seq.next() {
    //         let mut length = 0;
    //         while Some(current_base) == seq.peek().cloned() {
    //             seq.next();
    //             length += 1;
    //         }
    //         let index = match current_base {
    //             b'A' => 0,
    //             b'C' => 1,
    //             b'G' => 2,
    //             b'T' => 3,
    //             _ => panic!(),
    //         };
    //         if homo_polymers[index].len() <= length {
    //             let diff = length - homo_polymers[index].len() + 1;
    //             homo_polymers[index].extend(std::iter::repeat(0).take(diff));
    //         }
    //         homo_polymers[index][length] += 1;
    //     }
    // }
    // println!("Base\tLen\tCount");
    // for (base, counts) in homo_polymers.iter().enumerate() {
    //     for (len, count) in counts.iter().enumerate() {
    //         println!("{}\t{}\t{}", b"ACGT"[base] as char, len, count);
    //     }
    // }
    // use haplotyper::assemble::*;
    // let config = AssembleConfig::new(4, 1000, false, true, 6);
    // let min_count = 6;
    // let max_tig_length = 20;
    // ds.zip_up_suspicious_haplotig(&config, min_count, max_tig_length);
    // let target: u64 = args[2].parse().unwrap();
    // let reads: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter(|r| r.nodes.iter().any(|n| n.unit == target))
    // .collect();
    // let mut graph = ditch_graph::DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    // graph.remove_lightweight_edges(2, true);
    // let cov = ds.coverage.unwrap_or_else(|| panic!("Need coverage!"));
    // let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    // graph.clean_up_graph_for_assemble(cov, &lens, &reads, &config);
    // {
    //     graph.assign_copy_number(cov, &lens);
    //     graph.remove_zero_copy_elements(&lens, 0.2);
    //     graph.assign_copy_number(cov, &lens);
    //     graph.remove_zero_copy_elements(&lens, 0.5);
    //     graph.assign_copy_number(cov, &lens);
    //     graph.zip_up_overclustering();
    //     // From good Likelihood ratio focus, to weaker ones.
    //     for llr in (4..16).filter(|x| x % 4 == 0).rev() {
    //         graph.resolve_repeats(&reads, &config, llr as f64);
    //     }
    //     graph.zip_up_overclustering();
    //     graph.assign_copy_number(cov, &lens);
    //     graph.remove_zero_copy_path(0.3);
    // }
    // graph.remove_lightweight_edges(1, true);
    // for node in graph.nodes().filter(|n| n.node.0 == target) {
    //     println!("{}", node);
    // }
    // let min_count = 6;
    // let max_tig_length = 20;
    // ds.zip_up_suspicious_haplotig(&config, min_count, max_tig_length);
    // let config = AssembleConfig::new(4, 1000, false, false, 6);
    // println!("{}", ds.assemble(&config));

    // for read in ds.encoded_reads.iter() {
    //     for (node, edge) in read.nodes.iter().zip(read.edges.iter()) {
    //         assert_eq!(node.unit, edge.from);
    //     }
    // }
    // let is_ler: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.name.contains("ler")))
    //     .collect();
    // let target = (1713, 3151);
    // let ref_unit = ds.selected_chunks.iter().find(|n| n.id == 3302).unwrap();
    // let mut labels: Vec<_> = vec![];
    // for read in ds.encoded_reads.iter() {
    //     let is_ler = is_ler[&read.id];
    //     for edge in read.edges.iter() {
    //         let mut label = if (edge.from, edge.to) == target {
    //             edge.label.as_bytes().to_vec()
    //         } else if (edge.to, edge.from) == target {
    //             bio_utils::revcmp(edge.label.as_bytes())
    //         } else {
    //             continue;
    //         };
    //         label.iter_mut().for_each(u8::make_ascii_uppercase);
    //         let (_, ops) =
    //             kiley::bialignment::global_banded(&label, ref_unit.seq(), 2, -7, -4, -1, 100);
    //         let indel_iter = ops.iter().map(|op| match op {
    //             kiley::bialignment::Op::Mat => -1,
    //             _ => 1,
    //         });
    //         let dist = ops
    //             .iter()
    //             .filter(|&&op| !matches!(op, kiley::bialignment::Op::Mat))
    //             .count();
    //         use haplotyper::encode::max_region;
    //         let max_indel = max_region(indel_iter);
    //         println!("======={},{},{}========", is_ler, max_indel, dist);
    //         let (xr, ar, yr) = kiley::bialignment::recover(&label, ref_unit.seq(), &ops);
    //         for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //             println!("{}", String::from_utf8_lossy(xr));
    //             println!("{}", String::from_utf8_lossy(ar));
    //             println!("{}\n", String::from_utf8_lossy(yr));
    //         }
    //         labels.push(label);
    //     }
    // }
    // let med_idx = labels.len() / 2;
    // labels.select_nth_unstable_by_key(med_idx, |x| x.len());
    // labels.swap(0, med_idx);
    // let mut lengths: Vec<_> = labels.iter().map(|x| x.len()).collect();
    // let (_, median, _) = lengths.select_nth_unstable(labels.len() / 2);
    // let (upper, lower) = (2 * *median, (*median / 2).max(200));
    // let labels: Vec<&[u8]> = labels
    //     .iter()
    //     .filter(|ls| lower < ls.len() && ls.len() < upper)
    //     .map(|x| x.as_slice())
    //     .collect();
    // let rough_contig = kiley::ternary_consensus_by_chunk(&labels, 100);
    // let rough_contig =
    //     kiley::bialignment::polish_until_converge_banded(&rough_contig, &labels, 100);
    // let (_, ops) =
    //     kiley::bialignment::global_banded(&rough_contig, ref_unit.seq(), 2, -7, -4, -1, 100);
    // let (xr, ar, yr) = kiley::bialignment::recover(&rough_contig, ref_unit.seq(), &ops);
    // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //     println!("{}", String::from_utf8_lossy(xr));
    //     println!("{}", String::from_utf8_lossy(ar));
    //     println!("{}\n", String::from_utf8_lossy(yr));
    // }
    // let id2seq: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    // let mut edge_encoding_patterns: HashMap<_, _> = HashMap::new();
    // edge_encoding_patterns.insert(((1005, 0, false), (2362, 0, false)), vec![(433, 3471)]);
    // let seq = ds
    //     .selected_chunks
    //     .iter()
    //     .find(|u| u.id == 3417)
    //     .unwrap()
    //     .seq();
    // let mut edge_consensus = HashMap::new();
    // edge_consensus.insert(((1005, 0, false), (2362, 0, false)), seq.to_vec());
    // for read in ds.encoded_reads.iter() {
    //     let seq = &id2seq[&read.id];
    //     haplotyper::dense_encoding::fill_edges_by_new_units(
    //         read,
    //         seq,
    //         &edge_encoding_patterns,
    //         &edge_consensus,
    //     );
    // }
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|r| r.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    // use haplotyper::assemble::*;
    // let config = AssembleConfig::new(4, 1000, false, false, 6);
    // ds.squish_small_contig(&config, 15);
    // ds.zip_up_suspicious_haplotig(&config, 6, 25);
    // println!("{}", ds.assemble(&config));
    // let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
    // let sim_thr = 0.35;
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(1)
    //     .build_global()
    //     .unwrap();
    // let _ds =
    //     haplotyper::encode::deletion_fill::correct_unit_deletion(ds, &mut failed_trials, sim_thr);
    // let mut failed_trials: Vec<_> = vec![vec![]; ds.encoded_reads.len()];
    // ds = haplotyper::encode::deletion_fill::correct_unit_deletion(ds, &mut failed_trials, 0.35);

    // Local clustering.
    // let cov = ds.coverage.unwrap();
    // // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // let target: u64 = args[2].parse().unwrap();
    // let (seqs, ids): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| {
    //         //let is_hapa = id2desc[&r.id].contains("255v2");
    //         // let is_hapa = id2name[&r.id].contains("ler");
    //         let is_hapa = id2name[&r.id].contains("hapA");
    //         r.nodes
    //             .iter()
    //             .filter(|n| n.unit == target)
    //             .map(|n| (n.seq(), is_hapa as u8))
    //             .collect::<Vec<_>>()
    //     })
    //     .unzip();
    // let cl = ds
    //     .selected_chunks
    //     .iter()
    //     .find(|n| n.id == target)
    //     .unwrap()
    //     .cluster_num as u8;
    // let mut config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, cl, cov);
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(3424);
    // let (asn, _, score) =
    //     haplotyper::local_clustering::kmeans::clustering(&seqs, &mut rng, &mut config).unwrap();
    // eprintln!("{}", score);
    // eprintln!("{:?}\n{:?}", ids, asn);
    // let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    // for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
    //     let ref_unit = units[&node.unit];
    //     let indel_iter = node.cigar.iter().map(|&op| match op {
    //         Op::Del(l) | Op::Ins(l) => l as i32,
    //         Op::Match(l) => -(l as i32),
    //     });
    //     use haplotyper::encode::max_region;
    //     let max_indel = max_region(indel_iter);
    //     let (mut npos, mut rpos) = (0, 0);
    //     let (nodeseq, refseq) = (node.seq(), ref_unit.seq());
    //     let mut mism = 0;
    //     for op in node.cigar.iter() {
    //         match *op {
    //             Op::Match(l) => {
    //                 mism += nodeseq
    //                     .iter()
    //                     .skip(npos)
    //                     .take(l)
    //                     .zip(refseq.iter().skip(rpos).take(l))
    //                     .filter(|(x, y)| x != y)
    //                     .count();
    //                 rpos += l;
    //                 npos += l;
    //             }
    //             Op::Del(l) => rpos += l,
    //             Op::Ins(l) => npos += l,
    //         }
    //     }
    //     println!("{}\t{}\t{}\t{}", node.unit, node.cluster, max_indel, mism);
    // }
    // let seq: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    // let mut count: HashMap<_, u32> = units.keys().map(|&x| (x, 0)).collect();
    // for read in ds.encoded_reads.iter_mut() {
    //     read.nodes.retain(|node| {
    //         let cigar = &node.cigar;
    //         let indel_iter = cigar.iter().map(|&op| match op {
    //             Op::Del(l) | Op::Ins(l) => l as i32,
    //             Op::Match(l) => -(l as i32),
    //         });
    //         use haplotyper::encode::max_region;
    //         let max_indel = max_region(indel_iter);
    //         if 50 < max_indel {
    //             // let ref_unit = &units[&node.unit];
    //             // let (xr, ar, yr) = node.recover(ref_unit);
    //             // eprintln!("==================");
    //             // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //             //     eprintln!("{}", String::from_utf8_lossy(xr));
    //             //     eprintln!("{}", String::from_utf8_lossy(ar));
    //             //     eprintln!("{}\n", String::from_utf8_lossy(yr));
    //             // }
    //             let _ = count.entry(node.unit).and_modify(|x| *x += 1);
    //             false
    //         } else {
    //             true
    //         }
    //     });
    //     let mut nodes = vec![];
    //     nodes.append(&mut read.nodes);
    //     *read = haplotyper::encode::nodes_to_encoded_read(read.id, nodes, &seq[&read.id]).unwrap();
    // }
    // for (key, count) in count.iter().filter(|&(_, &c)| c > 0) {
    //     eprintln!("REMOVED\t{}\t{}", key, count);
    // }
    // print!("{}", serde_json::ser::to_string(&ds).unwrap());
    // use haplotyper::assemble::*;
    // let config = AssembleConfig::new(3, 1000, false, false, 4);
    // ds.squish_small_contig(&config, 30);
    // println!("{}", ds.assemble(&config));
    // em algorithm.
    // let target: u64 = args[2].parse().unwrap();
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // let repeat_num = 10;
    // let to_regularize = true;
    // let coverage_thr = 5;
    // let ref_unit = ds.selected_chunks.iter().find(|n| n.id == target).unwrap();
    // let unit_id = ref_unit.id;
    // let reads: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
    //     .collect();
    // let k = ref_unit.cluster_num;
    // let (new_clustering, score, _) = (1..=k)
    //     .flat_map(|k| std::iter::repeat(k).take(repeat_num))
    //     .enumerate()
    //     .map(|(i, k)| {
    //         use haplotyper::em_correction::*;
    //         let seed = unit_id * (i * k) as u64;
    //         let config = Config::new(repeat_num, seed, k, unit_id, to_regularize, coverage_thr);
    //         em_clustering(&reads, &config)
    //     })
    //     .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
    //     .unwrap();
    // log::debug!("SCORE\t{:.3}", score);
    // for &(id, position, cl) in new_clustering.iter() {
    //     let answer = id2desc[&id].contains("252v");
    //     let read = reads.iter().find(|r| r.id == id).unwrap();
    //     let prev = read.nodes[position].cluster;
    //     assert_eq!(read.nodes[position].unit, target);
    //     eprintln!("RESULT\t{}\t{}\t{}", answer, prev, cl);
    // }
    // use std::collections::HashSet;
    // let target = (848, 335);
    // let reads: HashSet<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter_map(|r| {
    //         r.edges
    //             .iter()
    //             .any(|e| ((e.from, e.to) == target) || ((e.to, e.from) == target))
    //             .then(|| r.id)
    //     })
    //     .collect();
    // haplotyper::encode::deletion_fill::correct_unit_deletion_selected(&mut ds, &reads, 0.25);
    // let ref_unit: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.seq())).collect();
    // let start = std::time::Instant::now();
    // for node in ds
    //     .encoded_reads
    //     .iter_mut()
    //     .flat_map(|r| r.nodes.iter_mut())
    //     .take(1000)
    // {
    //     let (_, ops) =
    //         kiley::bialignment::global_banded(node.seq(), ref_unit[&node.unit], 2, -5, -6, -1, 100);
    //     node.cigar = haplotyper::encode::compress_kiley_ops(&ops);
    // }
    // let end = std::time::Instant::now();
    // println!("{:?}", (end - start));
    // let nodes: usize = ds.encoded_reads.iter().map(|r| r.nodes.len()).sum();
    // println!("{}", (end - start).as_millis() * nodes as u128 / 1_000_000);
    // let mut counts: HashMap<_, u32> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     for node in read.nodes.iter().filter(|n| n.unit == 1077) {
    //         for mode in read.nodes.iter().filter(|n| n.unit == 1075) {
    //             *counts.entry((node.cluster, mode.cluster)).or_default() += 1;
    //         }
    //     }
    // }
    // println!("{:?}", counts);
    // let unit: u64 = args[2].parse().unwrap();
    // let (ids, nodes): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| {
    //         r.nodes
    //             .iter()
    //             .filter(|n| n.unit == unit)
    //             .map(|n| (r.id, n))
    //             .collect::<Vec<_>>()
    //     })
    //     .unzip();
    // let ref_unit = ds.selected_chunks.iter().find(|u| u.id == unit).unwrap();
    // let cl = ref_unit.cluster_num as u8;
    // let cov = 26f64;
    // log::debug!("CLNUM\t{}", cl);
    // // let cov = ds.coverage.unwrap();
    // let mut config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, cl, cov);
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(3434290824);
    // let start = std::time::Instant::now();
    // let seqs: Vec<_> = nodes.iter().map(|n| n.seq()).collect();
    // let (asn, _, score) =
    //     haplotyper::local_clustering::kmeans::clustering(&seqs, &mut rng, &mut config).unwrap();
    // let end = std::time::Instant::now();
    // use std::collections::HashMap;
    // log::debug!("SCORE\t{}\t{}", score, (end - start).as_secs());
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // for (id, ((readid, asn), node)) in ids.iter().zip(asn).zip(nodes).enumerate() {
    //     let ans = id2desc[readid].contains("000252v2") as usize;
    //     let indel = node.cigar.iter().map(|&op| match op {
    //         Op::Match(l) => -(l as i32),
    //         Op::Del(l) | Op::Ins(l) => l as i32,
    //     });
    //     let max_indel = haplotyper::encode::max_region(indel);
    //     let desc = id2desc[readid];
    //     log::debug!("ANSWER\t{}\t{}\t{}\t{}\t{}", id, ans, asn, max_indel, desc);
    // }
    Ok(())
}

// fn gen_cluster<R: Rng>(rng: &mut R, error_rate: f64, target: u64, max: u64) -> u64 {
//     match rng.gen_bool(error_rate) {
//         true => target,
//         false => (target + rng.gen_range(1..max)) % max,
//     }
// }

use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
// use std::collections::HashMap;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(340928);
    for _ in 0..10 {
        let xs = kiley::gen_seq::generate_seq(&mut rng, 2000);
        let ys = kiley::gen_seq::generate_seq(&mut rng, 2000);
        let (score, ops) = kiley::bialignment::global_banded(&xs, &ys, 2, -6, -5, -1, 20);
        let indel_mism = ops
            .iter()
            .map(|&op| 1 - 2 * (op == kiley::bialignment::Op::Mat) as i32);
        let max_indel = haplotyper::encode::max_region(indel_mism).max(0) as usize;
        println!("{}\t{}", score, max_indel);
    }
    // let (xr, ar, yr) = kiley::bialignment::recover(&xs, &ys, &ops);
    // for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
    //     eprintln!("{}", String::from_utf8_lossy(xr));
    //     eprintln!("{}", String::from_utf8_lossy(ar));
    //     eprintln!("{}\n", String::from_utf8_lossy(yr));
    // }
    // let args: Vec<_> = std::env::args().collect();
    // let mut ds: DataSet = std::fs::File::open(&args[1])
    //     .map(BufReader::new)
    //     .map(|x| serde_json::de::from_reader(x).unwrap())?;
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|r| r.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    // use haplotyper::assemble::*;
    // let config = AssembleConfig::new(4, 1000, false, true, 6);
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
    // let cov = ds.coverage.unwrap();
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // let target: u64 = args[2].parse().unwrap();
    // let (seqs, ids): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| {
    //         let is_hapa = id2desc[&r.id].contains("255v2");
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
    // let target = 1257;
    // let repeat_num = 10;
    // let to_regularize = true;
    // let coverage_thr = 3;
    // let ref_unit = ds.selected_chunks.iter().find(|n| n.id == target).unwrap();
    // let unit_id = ref_unit.id;
    // let reads: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
    //     .collect();
    // let k = ref_unit.cluster_num;

    // let (new_clustering, _, _) = (1..=k)
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
    //(/ (/ (/ (* 450 3000) 4) 30) 3600.0)
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
    // use haplotyper::em_correction::Context;
    // let cluster_num = 2;
    // let error_rate = 0.13;
    // let num = 20;
    // let span_prob = 0.1;
    // let seed = 2343909;
    // println!("lk\tnull\tspan\tseed\tnum\trand_index\tclusternum");
    // for num in 15..40 {
    //     for span_prob in (0..)
    //         .map(|i| 0.05 + 0.02 * i as f64)
    //         .take_while(|&x| x < 0.5)
    //     {
    //         for seed in 0..20 {
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    // let contexts: Vec<_> = (0..cluster_num)
    //     .flat_map(|i| {
    //         (0..num)
    //             .map(|id| {
    //                 let id = i * num + id;
    //                 let index = 0;
    //                 let unit = 1;
    //                 let cluster = 0; //gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                 let mut forward = vec![(4, 0)];
    //                 if rng.gen_bool(span_prob) {
    //                     let elm = gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                     forward.push((0, elm));
    //                 };
    //                 let mut backward = vec![(3, 0)];
    //                 if rng.gen_bool(span_prob) {
    //                     let elm = gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                     backward.push((2, elm));
    //                 };
    //                 Context::with_attrs(id, index, unit, cluster, forward, backward)
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect();
    // let (null_asn, null_lk, null_ofs) =
    //     haplotyper::em_correction::em_clustering_inner(&contexts, 1, &mut rng);
    // let (asns, lk, ofs) = haplotyper::em_correction::em_clustering_inner(&contexts, 2, &mut rng);
    // let answer: Vec<_> = (0..cluster_num)
    //     .flat_map(|i| vec![i as u8; num as usize])
    //     .collect();
    // let (pred, k) = if null_lk - null_ofs < lk - ofs {
    //     (asns, 2)
    // } else {
    //     (null_asn, 1)
    // };
    // let pred: Vec<_> = pred.iter().map(|x| x.2 as u8).collect();
    // let randindex = haplotyper::local_clustering::rand_index(&answer, &pred);
    // println!("NULL:{},k=2:{}", null_ofs, ofs);
    // println!(
    //     "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
    //     lk, null_lk, span_prob, seed, num, randindex, k
    // );
    //         }
    //     }
    // }
    Ok(())
}

// fn gen_cluster<R: Rng>(rng: &mut R, error_rate: f64, target: u64, max: u64) -> u64 {
//     match rng.gen_bool(error_rate) {
//         true => target,
//         false => (target + rng.gen_range(1..max)) % max,
//     }
// }

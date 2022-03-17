#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
use haplotyper::local_clustering::{local_clustering_selected, LocalClustering};
use kiley::recover;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for read in ds.encoded_reads.iter_mut() {
        let len = read.recover_raw_read().len();
        assert_eq!(len, read.original_length);
    }
    for read in ds.encoded_reads.iter_mut() {
        let idx = read.nodes.len() / 2;
        read.remove(idx);
    }
    // let (read_error_rate, unit_error_rate, sd) =
    //     haplotyper::encode::deletion_fill::estimate_error_rate(&ds, 0.15);
    // let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    // eprintln!("{}", sd);

    // use haplotyper::encode::nodes_to_encoded_read;
    // for read in ds.encoded_reads.iter_mut() {
    //     let mut seq = read.recover_raw_read();
    //     seq.iter_mut().for_each(u8::make_ascii_uppercase);
    //     let id = read.id;
    //     read.nodes.retain(|node| {
    //         let (_, aln, _) = node.recover(units[&node.unit]);
    //         let diff = aln.iter().filter(|&&x| x != b'|').count();
    //         let error_rate = diff as f64 / aln.len() as f64;
    //         let expected =
    //             read_error_rate[id as usize] + unit_error_rate[node.unit as usize].max(0f64);
    //         let threshold = expected + 10f64 * sd;
    //         error_rate < threshold
    //     });
    //     let mut nodes = Vec::with_capacity(read.nodes.len());
    //     nodes.append(&mut read.nodes);
    //     *read = nodes_to_encoded_read(read.id, nodes, &seq).unwrap();
    // }
    // ds.sanity_check();
    // println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // env_logger::init();
    // let args: Vec<_> = std::env::args().collect();
    // let mut ds: DataSet = std::fs::File::open(&args[1])
    //     .map(BufReader::new)
    //     .map(|r| serde_json::de::from_reader(r).unwrap())?;
    // use haplotyper::encode::deletion_fill;
    // deletion_fill::correct_unit_deletion(&mut ds, 0.15);
    // let dump = serde_json::ser::to_string(&ds).unwrap();
    // println!("{}", dump);
    // let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    // use rayon::prelude::*;
    // let error_rates: Vec<_> = ds
    //     .encoded_reads
    //     .par_iter()
    //     .flat_map(|r| {
    //         r.nodes
    //             .par_iter()
    //             .map(|n| (r.id, n.unit, error(n, ref_units[&n.unit])))
    //             .collect::<Vec<_>>()
    //     })
    //     .collect();
    // println!("readid\tunitid\tmism\tins\tdel");
    // for (rid, unit, (mism, ins, del)) in error_rates {
    //     println!("{}\t{}\t{}\t{}\t{}", rid, unit, mism, ins, del);
    // }
    // let nodes: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| r.nodes.iter())
    //     .filter(|n| n.unit == 141)
    //     .collect();
    // let ref_unit = ds
    //     .selected_chunks
    //     .iter()
    //     .find(|c| c.id == 141)
    //     .unwrap()
    //     .seq();
    // let (seqs, mut ops): (Vec<_>, Vec<_>) = nodes
    //     .iter()
    //     .map(|node| {
    //         let ops = kiley::bialignment::global(ref_unit, node.seq(), 2, -6, -5, -2).1;
    //         (node.seq(), ops)
    //     })
    //     .unzip();
    // let cons =
    //     kiley::bialignment::guided::polish_until_converge_with(ref_unit, &seqs, &mut ops, 60);
    // for (seq, op) in seqs.iter().zip(ops.iter()) {
    //     let (refr, aln, query) = kiley::recover(&cons, seq, op);
    //     for ((refr, aln), query) in refr.chunks(200).zip(aln.chunks(200)).zip(query.chunks(200)) {
    //         println!("{}", String::from_utf8_lossy(refr));
    //         println!("{}", String::from_utf8_lossy(aln));
    //         println!("{}\n", String::from_utf8_lossy(query));
    //     }
    // }
    // let id2desc: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.name.clone()))
    //     .collect();
    // let mut nodes: Vec<_> = vec![];
    // let mut answer = vec![];
    // // let target = (210, 0);
    // for read in ds.encoded_reads.iter() {
    //     let is_hap1 = id2desc[&read.id].contains("000251v2") as usize;
    //     // for node in read.nodes.iter().filter(|n| (n.unit, n.cluster) == target) {
    //     for node in read.nodes.iter().filter(|n| n.unit == 210) {
    //         nodes.push(node.seq());
    //         answer.push(is_hap1);
    //     }
    // }
    // let mut rng: rand_xoshiro::Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(4329890);
    // use haplotyper::local_clustering::kmeans::ClusteringConfig;
    // let config = ClusteringConfig::new(100, 4, 25f64, definitions::ReadType::CLR);
    // let (preds, gains, _, _) =
    //     haplotyper::local_clustering::kmeans::clustering(&nodes, &mut rng, &config).unwrap();
    // for ((p, a), g) in preds.iter().zip(answer.iter()).zip(gains.iter()) {
    //     println!("{}\t{}\t{}", a, p, vec2str(g));
    // }
    Ok(())
}

// fn vec2str(xs: &[f64]) -> String {
//     let xs: Vec<_> = xs.iter().map(|&x| format!("{:6.1}", x)).collect();
//     xs.join(",")
// }

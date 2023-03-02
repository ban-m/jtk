use definitions::*;
const BIAS_THR: f64 = 0.2;
use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let mut chunk_pairs: HashMap<_, usize> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let nodes = read.nodes.iter().enumerate();
        for (i, n1) in nodes.filter(|n| n.1.is_biased(BIAS_THR)) {
            let n2s = read.nodes.iter().skip(i + 1);
            for n2 in n2s.filter(|n| n.is_biased(BIAS_THR)) {
                let key = (n1.chunk.min(n2.chunk), n1.chunk.max(n2.chunk));
                *chunk_pairs.entry(key).or_default() += 1;
            }
        }
    }
    let chunks: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|n| (n.id, n.cluster_num))
        .collect();
    chunk_pairs.retain(|_, val| 8 < *val);
    chunk_pairs.retain(|(u1, u2), _| 1 < chunks[u1] && 1 < chunks[u2]);
    let mut max_correl: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, 0f64)).collect();
    for &(u1, u2) in chunk_pairs.keys() {
        let (cl1, cl2) = (chunks[&u1], chunks[&u2]);
        let (rel, _) = check_correl(&ds, (u1, cl1), (u2, cl2));
        max_correl.entry(u1).and_modify(|x| *x = x.max(rel));
        max_correl.entry(u2).and_modify(|x| *x = x.max(rel));
    }
    let sum: f64 = max_correl.values().sum();
    let ave = sum / max_correl.len() as f64;
    println!("{},{:.3},{:.3}", max_correl.len(), sum, ave);
    Ok(())
}

fn check_correl(
    ds: &DataSet,
    (chunk1, cl1): (u64, usize),
    (chunk2, cl2): (u64, usize),
) -> (f64, usize) {
    let (mut c1, mut c2) = (vec![], vec![]);
    for read in ds.encoded_reads.iter() {
        let node1 = read
            .nodes
            .iter()
            .filter(|n| n.chunk == chunk1 && n.is_biased(BIAS_THR))
            .map(|n| n.cluster as usize)
            .min();
        let node2 = read
            .nodes
            .iter()
            .filter(|n| n.chunk == chunk2 && n.is_biased(BIAS_THR))
            .map(|n| n.cluster as usize)
            .min();
        if let (Some(n1), Some(n2)) = (node1, node2) {
            c1.push(n1);
            c2.push(n2);
        }
    }
    if c1.is_empty() {
        return (0f64, c1.len());
    }
    let c1_is_same = c1.iter().all(|&x| x == c1[0]);
    let c2_is_same = c2.iter().all(|&x| x == c2[0]);
    let rel_value = match (c1_is_same && c2_is_same, cl1 == 1 && cl2 == 1) {
        (true, true) => 0f64,
        (true, false) => 1f64,
        (false, _) => haplotyper::misc::adjusted_rand_index(&c1, &c2),
    };
    if rel_value.is_nan() {
        return (0f64, c1.len());
    }
    (rel_value, c1.len())
}

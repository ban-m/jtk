use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let unit = ds.selected_chunks.iter().find(|n| n.id == 1707).unwrap();
    let dists: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|r| r.nodes.iter().filter(|n| n.unit == 1707))
        .map(|node| edlib_sys::global_dist(unit.seq(), node.seq()))
        .collect();
    for d in dists.iter() {
        println!("{}", d);
    }
    let (sum, sumsq) = dists
        .iter()
        .fold((0, 0), |(sum, sumsq), x| (sum + x, sumsq + x * x));
    let mean = sum as f64 / dists.len() as f64;
    let var = sumsq as f64 / dists.len() as f64 - mean * mean;
    println!("{}\t{}", mean, var.sqrt());
    // let mut edge_count: HashMap<_, usize> = HashMap::new();
    // let mut edge_lens: HashMap<_, i64> = HashMap::new();
    // for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
    //     let len = if edge.label.is_empty() {
    //         edge.offset
    //     } else {
    //         edge.label().len() as i64
    //     };
    //     let edge = (edge.from.min(edge.to), edge.from.max(edge.to));
    //     *edge_count.entry(edge).or_default() += 1;
    //     *edge_lens.entry(edge).or_default() += len;
    // }
    // println!("From\tTo\tCount");
    // for (&(from, to), &value) in edge_count.iter() {
    //     let len = edge_lens[&(from, to)] / value as i64;
    //     println!("{}\t{}\t{}\t{}", from, to, value, len);
    // }
    Ok(())
}

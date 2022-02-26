#![allow(unused_imports)]
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for r in ds.encoded_reads.iter() {
        println!("{}\t{}", r.id, r.original_length);
    }
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // let mut nodes: HashMap<u64, Vec<(&Node, bool)>> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     let ans = match IS_MOCK {
    //         true => id2desc[&read.id].contains("hapA"),
    //         false => id2desc[&read.id].contains("000251v2"),
    //     };
    //     for node in read.nodes.iter() {
    //         nodes.entry(node.unit).or_default().push((node, ans));
    //     }
    // }
    // let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    // println!("Unit\tCoverage\tRndIndex\tCopyNum\tLK");
    // for (id, nodes) in nodes {
    //     let chunk = chunks[&id];
    //     let answer: Vec<_> = nodes.iter().map(|x| x.1 as u8).collect();
    //     let preds: Vec<_> = nodes.iter().map(|x| x.0.cluster as u8).collect();
    //     let rand = haplotyper::local_clustering::rand_index(&answer, &preds);
    //     println!(
    //         "{}\t{}\t{}\t{}\t{}",
    //         id,
    //         nodes.len(),
    //         rand,
    //         chunk.copy_num,
    //         chunk.score,
    //     );
    // }
    Ok(())
}

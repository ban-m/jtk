const IS_MOCK: bool = false;
use definitions::*;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    let reads: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    for read in ds.encoded_reads.iter().filter(|r| reads.contains(&r.id)) {
        let dump: Vec<_> = read
            .nodes
            .iter()
            .map(|n| (n.unit, n.cluster, n.is_forward, n.position_from_start))
            .collect();
        let is_hap1 = match IS_MOCK {
            false => id2desc[&read.id].contains("251v2"),
            true => id2desc[&read.id].contains("hapA"),
        };
        println!("{}\t{}\t{:?}", read.id, is_hap1, dump);
    }
    //let units: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    // let range = 10;
    // for read in ds.encoded_reads.iter() {
    //     let mut dumps = vec![format!("{:<5}", 0); range];
    //     for (idx, node) in read
    //         .nodes
    //         .iter()
    //         .enumerate()
    //         .filter(|(_, n)| units.contains(&n.unit))
    //     {
    // let (nodes, idx) = {
    //     let mut nodes: Vec<_> = read.nodes.iter().map(|n| n.unit).collect();
    //     match node.is_forward {
    //         true => (nodes, idx),
    //         false => {
    //             nodes.reverse();
    //             let idx = nodes.len() - idx - 1;
    //             (nodes, idx)
    //         }
    //     }
    // };
    // for i in (0..range).filter(|&i| range / 2 <= idx + i) {
    //     if let Some(unit) = nodes.get(idx + i - range / 2) {
    //         dumps[i] = format!("{:<5}", unit);
    //     }
    // }
    // let is_hap1 = match IS_MOCK {
    //     false => id2desc[&read.id].contains("251v2"),
    //     true => id2desc[&read.id].contains("hapA"),
    // };
    // println!("{}\t{}\t{}", is_hap1, node.cluster, dumps.join("\t"));
    // println!(
    //     "{}\t{}\t{}\t{}",
    //     read.id,
    //     is_hap1,
    //     node.cluster,
    //     dumps.join("\t")
    // );
    //     }
    // }
    Ok(())
}

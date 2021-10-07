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
        .map(|r| (r.id, r.desc.clone()))
        .collect();
    let units: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    for read in ds.encoded_reads.iter() {
        let mut dumps = vec![format!("{:<5}", 0); 7];
        for (idx, node) in read
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| units.contains(&n.unit))
        {
            let (nodes, idx) = {
                let mut nodes: Vec<_> = read.nodes.iter().map(|n| n.unit).collect();
                match node.is_forward {
                    true => (nodes, idx),
                    false => {
                        nodes.reverse();
                        let idx = nodes.len() - idx - 1;
                        (nodes, idx)
                    }
                }
            };
            for i in (0..7).filter(|&i| 3 <= idx + i) {
                if let Some(unit) = nodes.get(idx + i - 3) {
                    dumps[i] = format!("{:<5}", unit);
                }
            }
            let is_hap1 = id2desc[&read.id].contains("252v2");
            // println!("{}\t{}\t{}", is_hap1, node.cluster, dumps.join("\t"));
            println!(
                "{}\t{}\t{}\t{}",
                read.id,
                is_hap1,
                node.cluster,
                dumps.join("\t")
            );
        }
    }
    Ok(())
}

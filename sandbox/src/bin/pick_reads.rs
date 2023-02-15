use definitions::*;
use sandbox::IS_MOCK;
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

    let chunks: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    let range = 14;
    for read in ds.encoded_reads.iter() {
        let len = read.original_length;
        for (idx, node) in read
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| chunks.contains(&n.chunk))
        {
            let mut dumps = vec![format!("{:<5}", 0); range];
            let (nodes, idx) = {
                let mut nodes: Vec<_> = read.nodes.iter().map(|n| n.chunk).collect();
                match node.is_forward {
                    true => (nodes, idx),
                    false => {
                        nodes.reverse();
                        let idx = nodes.len() - idx - 1;
                        (nodes, idx)
                    }
                }
            };
            for i in (0..range).filter(|&i| range / 2 <= idx + i) {
                if let Some(chunk) = nodes.get(idx + i - range / 2) {
                    dumps[i] = format!("{:<5}", chunk);
                }
            }
            let is_hap1 = match IS_MOCK {
                false => id2desc[&read.id].contains("251v2"),
                true => id2desc[&read.id].contains("hapA"),
            };
            let dir = node.is_forward;
            println!(
                "{is_hap1}\t{dir}\t{}\t{len}\t{}\t{}",
                node.cluster,
                dumps.join("\t"),
                read.id,
            );
        }
    }
    Ok(())
}

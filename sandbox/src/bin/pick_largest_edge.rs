use definitions::*;
// use std::collections::HashMap;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let targets: Vec<_> = args[2..]
        .iter()
        .map(|x| {
            let (chunk, cluster) = x.rsplit_once('-').unwrap();
            let chunk: u64 = chunk.parse().unwrap();
            let cluster: u64 = cluster.parse().unwrap();
            (chunk, cluster)
        })
        .collect();
    for &target in targets.iter() {
        println!("=========={target:?}==========");
        for read in ds.encoded_reads.iter() {
            for (i, node) in read.nodes.iter().enumerate() {
                if (node.chunk, node.cluster) == target {
                    if 0 < i {
                        let edge = &read.edges[i - 1];
                        let seq = std::str::from_utf8(edge.label()).unwrap();
                        let to = (read.nodes[i - 1].chunk, read.nodes[i - 1].cluster);
                        println!("{to:?}\t{}\t{seq}", seq.len());
                    }
                    if i + 1 < read.nodes.len() {
                        let edge = &read.edges[i];
                        let seq = std::str::from_utf8(edge.label()).unwrap();
                        let to = &read.nodes[i + 1];
                        let to = (to.chunk, to.cluster);
                        println!("{to:?}\t{}\t{seq}", seq.len());
                    }
                }
            }
        }
    }

    Ok(())
}

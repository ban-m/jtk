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
            let (unit, cluster) = x.rsplit_once('-').unwrap();
            let unit: u64 = unit.parse().unwrap();
            let cluster: u64 = cluster.parse().unwrap();
            (unit, cluster)
        })
        .collect();
    for &target in targets.iter() {
        println!("=========={target:?}==========");
        for read in ds.encoded_reads.iter() {
            for (i, node) in read.nodes.iter().enumerate() {
                if (node.unit, node.cluster) == target {
                    if 0 < i {
                        let edge = &read.edges[i - 1];
                        let seq = std::str::from_utf8(edge.label()).unwrap();
                        let to = (read.nodes[i - 1].unit, read.nodes[i - 1].cluster);
                        println!("{to:?}\t{}\t{seq}", seq.len());
                    }
                    if i + 1 < read.nodes.len() {
                        let edge = &read.edges[i];
                        let seq = std::str::from_utf8(edge.label()).unwrap();
                        let to = &read.nodes[i + 1];
                        let to = (to.unit, to.cluster);
                        println!("{to:?}\t{}\t{seq}", seq.len());
                    }
                }
            }
        }
    }
    // let mut edge_count: HashMap<_, Vec<_>> = HashMap::new();
    // for read in ds.encoded_reads.iter().filter(|r| !r.nodes.is_empty()) {
    //     assert_eq!(read.nodes.len(), read.edges.len() + 1);
    //     for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
    //         let label = edge.label();
    //         let (edge, _is_forward) = normalize_edge(w);
    //         edge_count.entry(edge).or_default().push(label);
    //     }
    // }
    // edge_count.retain(|_, edges| edges.len() < 3);
    // let (max, argmax) = edge_count
    //     .iter()
    //     .max_by_key(|(_, edges)| edges.iter().map(|x| x.len()).max().unwrap())
    //     .unwrap();
    // eprintln!("{max:?}");
    // let seq = argmax.iter().max_by_key(|x| x.len()).unwrap();
    // eprintln!("{}", std::str::from_utf8(seq).unwrap());
    Ok(())
}

// type NormedEdge = ((u64, bool), (u64, bool));
// fn normalize_edge(w: &[Node]) -> (NormedEdge, bool) {
//     let forward = ((w[0].unit, w[0].is_forward), (w[1].unit, w[1].is_forward));
//     let reverse = ((w[1].unit, !w[1].is_forward), (w[0].unit, !w[0].is_forward));
//     if forward <= reverse {
//         (forward, true)
//     } else {
//         (reverse, false)
//     }
// }

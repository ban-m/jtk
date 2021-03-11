use definitions::*;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // let units: HashSet<u64> = std::fs::File::open(&args[2])
    //     .map(BufReader::new)?
    //     .lines()
    //     .filter_map(|x| x.ok())
    //     .filter_map(|x| {
    //         // let mut x = x.split('-');
    //         // let unit: u64 = x.next()?.parse().ok()?;
    //         // let cluster: u64 = x.next()?.parse().ok()?;
    //         // Some((unit, cluster))
    //         x.parse().ok()
    //     })
    let units: HashSet<u64> = vec![1635].into_iter().collect();
    use std::collections::HashMap;
    let id2name: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.name.to_string()))
        .collect();
    for read in ds.encoded_reads.iter() {
        if read
            .nodes
            .iter()
            // .any(|n| units.contains(&(n.unit, n.cluster)))
            .any(|n| units.contains(&n.unit))
        {
            let mut line = vec![];
            for (w, e) in read.nodes.windows(2).zip(read.edges.iter()) {
                line.push(format!("{}({})", w[0].unit, e.offset));
            }
            println!("{}\t{}", id2name[&read.id], line.join(""));
        }
    }
    Ok(())
}

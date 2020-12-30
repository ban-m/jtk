use definitions::*;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let units: HashSet<(u64, u64)> = std::fs::File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|x| {
            let mut x = x.split('-');
            let unit: u64 = x.next()?.parse().ok()?;
            let cluster: u64 = x.next()?.parse().ok()?;
            Some((unit, cluster))
        })
        .collect();
    for read in ds.encoded_reads.iter() {
        if read
            .nodes
            .iter()
            .any(|n| units.contains(&(n.unit, n.cluster)))
        {
            let mut line = vec![];
            for (w, e) in read.nodes.windows(2).zip(read.edges.iter()) {
                line.push(format!("{}({}){}", w[0].unit, e.offset, w[1].unit));
            }
            println!("{}", line.join("-"));
        }
    }
    Ok(())
}

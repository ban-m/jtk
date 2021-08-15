use definitions::*;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    println!("ID1\tID2\tK1\tK2\tCount\tPValue");
    let mut pair_units: HashMap<_, HashMap<_, u32>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for (i, n) in read.nodes.iter().enumerate() {
            for m in read.nodes.iter().skip(i + 1) {
                let (n, m) = (n.unit.min(m.unit), n.unit.max(m.unit));
                *(pair_units.entry(n).or_default().entry(m).or_default()) += 1;
            }
        }
    }
    for (id, connections) in pair_units {
        for (&jd, count) in connections.iter().filter(|&(_, &count)| 10 < count) {
            let mut pairs = vec![];
            for read in ds.encoded_reads.iter() {
                let read = read.nodes.iter().enumerate();
                for (idx, node) in read.clone().filter(|(_, n)| n.unit == id) {
                    let cluster = node.cluster;
                    for (_, n) in read.clone().filter(|&(i, n)| i != idx && n.unit == jd) {
                        pairs.push((cluster, n.cluster));
                    }
                }
            }
            let p_value = haplotyper::unit_correlation::calc_p_value(&pairs, 4324);
            let (cluster1, cluster2): (HashSet<_>, HashSet<_>) = pairs.iter().copied().unzip();
            let (k1, k2) = (cluster1.len(), cluster2.len());
            println!("{}\t{}\t{}\t{}\t{}\t{}", id, jd, k1, k2, count, p_value);
            println!("{}\t{}\t{}\t{}\t{}\t{}", jd, id, k2, k1, count, p_value);
        }
    }
    Ok(())
}

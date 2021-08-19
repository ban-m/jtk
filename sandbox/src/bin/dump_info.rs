use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use std::collections::HashMap;
    let mut counts: HashMap<(u64, u64), u32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *counts.entry((node.unit, node.cluster)).or_default() += 1;
        }
    }
    println!("UNIT\tunit\tcluster\tcount");
    for ((unit, cluster), counts) in counts.iter() {
        println!("UNIT\t{}\t{}\t{}", unit, cluster, counts);
    }
    use std::cmp::Ordering;
    let edge = |(u1, c1): (u64, u64), (u2, c2): (u64, u64)| match u1.cmp(&u2) {
        Ordering::Greater => ((u2, c2), (u1, c1)),
        Ordering::Less => ((u1, c1), (u2, c2)),
        Ordering::Equal => match c1.cmp(&c2) {
            Ordering::Greater => ((u2, c2), (u1, c1)),
            Ordering::Less => ((u1, c1), (u2, c2)),
            _ => ((u1, c1), (u2, c2)),
        },
    };
    let mut counts: HashMap<((u64, u64), (u64, u64)), (u32, i64)> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for w in read.nodes.windows(2) {
            let edge_len = w[1].position_from_start as i64
                - (w[0].position_from_start + w[0].query_length()) as i64;
            let w1 = (w[0].unit, w[0].cluster);
            let w2 = (w[1].unit, w[1].cluster);
            let _ = *counts
                .entry(edge(w1, w2))
                .and_modify(|(count, len)| {
                    *count += 1;
                    *len += edge_len;
                })
                .or_insert((1, edge_len));
        }
    }
    println!("EDGE\tfrom\tfcluster\tto\ttcluster\tcount\tlen");
    for (((u1, c1), (u2, c2)), (c, len)) in counts {
        let len = len / c as i64;
        println!("EDGE\t{}\t{}\t{}\t{}\t{}\t{}", u1, c1, u2, c2, c, len);
    }
    Ok(())
}

use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let mut counts: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    for (key, val) in counts.iter() {
        println!("{}\t{}", key, val);
    }
    let mut counts: Vec<_> = counts.values().copied().collect();
    counts.sort();
    eprintln!("{}", counts[counts.len() / 2]);
    // use haplotyper::re_clustering::*;
    // let unit_id: u64 = args[2].parse().unwrap();
    // use haplotyper::local_clustering::kmeans;
    // use rand::SeedableRng;
    // let mut rng: rand_xoshiro::Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(unit_id * 23);
    // let units: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| r.nodes.iter().filter(|n| n.unit == unit_id))
    //     .collect();
    // let cluster_num = 12;
    // let seqs: Vec<_> = units.iter().map(|n| n.seq()).collect();
    // let coverage = ds.coverage.unwrap();
    // log::debug!("Coverage\t{}", coverage);
    // let mut config = kmeans::ClusteringConfig::new(100, cluster_num, coverage);
    // let (asn, _consensus) = kmeans::clustering(&seqs, &mut rng, &mut config).unwrap();
    // let ids: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|r| {
    //         let len = r.nodes.iter().filter(|n| n.unit == unit_id).count();
    //         vec![r.id; len]
    //     })
    //     .collect();
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // for (id, asn) in ids.iter().zip(asn) {
    //     let answer = id2desc[id].contains("000252v");
    //     log::debug!("RESULT\t{}\t{}\t{}", id, asn, answer);
    // }
    // let target = 1374;
    // for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
    //     if edge.from == target || edge.to == target {
    //         let (x, y) = if edge.from == target {
    //             (edge.from, edge.to)
    //         } else {
    //             (edge.to, edge.from)
    //         };
    //         println!("{}\t{}\t{}", x, y, edge.offset);
    //     }
    // }
    // let unit = ds.selected_chunks.iter().find(|u| u.id == 375).unwrap();
    // for read in ds.encoded_reads.iter() {
    //     for (idx, node) in read.nodes.iter().enumerate().filter(|(_, n)| n.unit == 375) {
    //         let (_, ax, _) = node.recover(&unit);
    //         let dist = ax.iter().filter(|&&x| x != b'|').count();
    //         let to = match node.is_forward {
    //             true => read.nodes.get(idx + 1).map(|x| x.unit).unwrap_or(0),
    //             false => read.nodes.get(idx - 1).map(|x| x.unit).unwrap_or(0),
    //         };
    //         println!("DIST\t{}\t{}\t{}", unit.id, to, dist);
    // for ((qx, ax), rx) in qx.chunks(200).zip(ax.chunks(200)).zip(rx.chunks(200)) {
    //     let qx = String::from_utf8_lossy(qx);
    //     let ax = String::from_utf8_lossy(ax);
    //     let rx = String::from_utf8_lossy(rx);
    //     println!("{}\n{}\n{}\n", qx, ax, rx);
    // }
    //     }
    // }
    Ok(())
}

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
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        let ref_chunk = chunks.get(&node.unit).unwrap();
        let (_, aln, _) = node.recover(ref_chunk);
        let dist = aln.iter().filter(|&&x| x != b'|').count() as f64 / aln.len() as f64;
        println!("{}\t{}", node.unit, dist);
    }
    // let mut counts: HashMap<(u64, u64), u32> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     for node in read.nodes.iter() {
    //         *counts.entry((node.unit, node.cluster)).or_default() += 1;
    //     }
    // }
    // let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
    // println!("UNIT\tunit\tcluster\tcount\tscore");
    // for ((unit, cluster), counts) in counts.iter() {
    //     let score = score[unit];
    //     println!("UNIT\t{}\t{}\t{}\t{:.4}", unit, cluster, counts, score);
    // }
    // let mut edge_counts: HashMap<((u64, u64), (u64, u64)), u32> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     for w in read.nodes.windows(2) {
    //         let from = (w[0].unit, w[0].cluster);
    //         let to = (w[1].unit, w[1].cluster);
    //         *edge_counts.entry((from.min(to), to.max(from))).or_default() += 1;
    //     }
    // }
    // println!("EDGE\tfrom\tfromcl\tto\ttocl\tcount");
    // for (((from, fromcl), (to, tocl)), count) in edge_counts {
    //     println!("EDGE\t{}\t{}\t{}\t{}\t{}", from, fromcl, to, tocl, count);
    // }
    Ok(())
}

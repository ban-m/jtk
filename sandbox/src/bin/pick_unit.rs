use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
// Synopsis: [JSON] [Unit to dump]+
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
    for unit in args[2..].iter().map(|x| -> u64 { x.parse().unwrap() }) {
        for read in ds.encoded_reads.iter() {
            // let is_hap1 = id2desc[&read.id].contains("251v2");
            let is_hap1 = id2desc[&read.id].contains("hapA");
            for node in read.nodes.iter().filter(|n| unit == n.unit) {
                let post: Vec<_> = node.posterior.iter().map(|p| format!("{:.3}", p)).collect();
                let (unit, cluster, post) = (node.unit, node.cluster, post.join("\t"));
                println!("{}\t{}\t{}\t{}", is_hap1, unit, cluster, post);
            }
        }
    }
    Ok(())
}

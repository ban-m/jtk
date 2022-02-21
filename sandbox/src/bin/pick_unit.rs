const IS_MOCK: bool = false;
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
        let ref_chunk = ds.selected_chunks.iter().find(|c| c.id == unit).unwrap();
        for read in ds.encoded_reads.iter() {
            let is_hap1 = match IS_MOCK {
                true => id2desc[&read.id].contains("hapA") as usize,
                false => id2desc[&read.id].contains("000251v2") as usize,
            };
            for (i, node) in read.nodes.iter().enumerate().filter(|n| unit == n.1.unit) {
                let (query, aln, refr) = node.recover(ref_chunk);
                let dist = aln.iter().filter(|&&x| x != b'|').count();
                let identity = 1f64 - dist as f64 / aln.len() as f64;
                let post: Vec<_> = node.posterior.iter().map(|p| format!("{:.3}", p)).collect();
                let (unit, cluster, post) = (node.unit, node.cluster, post.join("\t"));
                let len = read.nodes.len();
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                    is_hap1, i, len, unit, cluster, identity, post
                );
                println!("ALN\t{}", dist);
                for (query, refr) in query.chunks(200).zip(refr.chunks(200)) {
                    println!("ALN\t{}", String::from_utf8_lossy(query));
                    println!("ALN\t{}", String::from_utf8_lossy(refr));
                }
            }
        }
    }
    Ok(())
}

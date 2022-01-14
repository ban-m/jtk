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
    Ok(())
}

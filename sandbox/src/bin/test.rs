use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2name: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.name.clone()))
        .collect();
    println!("ID\tScore");
    for unit in ds.selected_chunks.iter() {
        let (pred, answer): (Vec<_>, Vec<_>) = ds
            .encoded_reads
            .iter()
            .flat_map(|r| {
                let answer = id2name[&r.id].starts_with("hapA");
                r.nodes
                    .iter()
                    .filter(|n| n.unit == unit.id)
                    .map(|n| (n.cluster as u8, answer as u8))
                    .collect::<Vec<_>>()
            })
            .unzip();
        let score = haplotyper::local_clustering::rand_index(&pred, &answer);
        println!("{}\t{}", unit.id, score);
    }
    Ok(())
}

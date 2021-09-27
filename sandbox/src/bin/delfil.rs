use definitions::*;
use std::collections::HashSet;
use std::io::*;
// Fill deletion test.
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    // let target = 520;
    // let mut neighbor: HashSet<_> = HashSet::new();
    // for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
    //     if edge.from == target {
    //         neighbor.insert(edge.to);
    //     }
    //     if edge.to == target {
    //         neighbor.insert(edge.from);
    //     }
    // }
    // let target_reads: HashSet<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter_map(|r| {
    //         r.nodes
    //             .iter()
    //             .any(|n| neighbor.contains(&n.unit))
    //             .then(|| r.id)
    //     })
    //     .collect();
    let mut failed = vec![vec![]; ds.encoded_reads.len()];
    ds = haplotyper::encode::deletion_fill::correct_unit_deletion(ds, &mut failed, 0.35);
    Ok(())
}

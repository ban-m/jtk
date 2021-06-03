use definitions::*;
use log::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    debug!("{:?}", end - start);
    // use std::collections::HashMap;
    // use haplotyper::GlobalClustering;
    // let config = haplotyper::GlobalClusteringConfig::new(3, 10, 1, -1, 1);
    // ds = ds.global_clustering_graph(&config);
    // let id2name: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|read| (read.id, &read.name))
    //     .collect();
    // for &Assignment { id, cluster } in ds.assignments.iter() {
    //     println!("{}\t{}\t{}", id, id2name[&id], cluster);
    // }
    // use haplotyper::ComponentPicking;
    // let config = haplotyper::ComponentPickingConfig::new(1);
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    // ds.encoded_reads
    //     .iter_mut()
    //     .for_each(|r| r.nodes.iter_mut().for_each(|n| n.cluster = 0));
    // let remove_unit = vec![376, 1648];
    // let remove_unit: Vec<_> = ds
    //     .selected_chunks
    //     .iter()
    //     .filter_map(|u| (4 <= u.cluster_num).then(|| u.id))
    //     .collect();
    // for read in ds.encoded_reads.iter_mut() {
    //     while let Some(idx) = read
    //         .nodes
    //         .iter()
    //         .position(|n| remove_unit.contains(&n.unit))
    //     {
    //         read.remove(idx);
    //     }
    // }
    let config = haplotyper::assemble::AssembleConfig::new(24, 100, false);
    use haplotyper::Assemble;
    let gfa = ds.assemble_as_gfa(&config);
    println!("{}", gfa);
    Ok(())
}

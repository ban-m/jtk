use definitions::*;
use log::*;
// use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    debug!("{:?}", end - start);
    let config = haplotyper::assemble::AssembleConfig::new(24, 100, false);
    use haplotyper::Assemble;
    ds.assignments = ds
        .raw_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    let gfa = ds.assemble_as_gfa(&config);
    println!("{}", gfa);
    Ok(())
    // use haplotyper::assemble::string_graph::*;
    // let config = AssembleConfig::new(2, 0.2);
    // let mut graph = StringGraph::new(&ds.encoded_reads, &config);
    // graph.node_count();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.transitive_edge_reduction();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.simple_path_reduction();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // let (nodes, _, group, summaries) = graph.assemble();
    // let nodes_length: HashMap<_, _> = nodes.iter().map(|n| (n.sid.clone(), n.slen)).collect();
    // let retained_reads: HashSet<_> = {
    //     let mut retained_reads = HashSet::new();
    //     for g in group {
    //         let g = g.set().unwrap();
    //         let is_short_contig = g.iter().all(|id| nodes_length.get(id).unwrap() < &100_000);
    //         if !is_short_contig || 1 < g.iter().count() {
    //             for contig in g.iter() {
    //                 let summary = summaries.iter().find(|summary| &summary.id == contig);
    //                 for read in summary.unwrap().reads.iter() {
    //                     retained_reads.insert(read);
    //                 }
    //             }
    //         }
    //     }
    //     retained_reads
    // };
    // ds.encoded_reads.retain(|r| retained_reads.contains(&r.id));
    // let config = AssembleConfig::new(2, 0.2);
    // let mut graph = StringGraph::from_dataset(&ds, &config);
    // graph.node_count();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.transitive_edge_reduction();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.simple_path_reduction();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.light_node_trim(5);
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // graph.simple_path_reduction();
    // assert!(graph.sanity_check());
    // debug!("{:?}", graph);
    // let gfa = graph.assemble_as_gfa();
    // println!("{}", gfa);
    // Ok(())
}

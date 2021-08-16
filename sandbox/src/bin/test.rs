use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use haplotyper::assemble::ditch_graph::DitchGraph;
    use haplotyper::assemble::AssembleConfig;
    let config = AssembleConfig::new(2, 100, false, true);
    {
        let reads: Vec<_> = ds.encoded_reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
        graph.remove_lightweight_edges(2);
        let squish = graph.squish_bubbles(15);
        ds.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .for_each(|n| match squish.get(&(n.unit, n.cluster)) {
                Some(res) => n.cluster = *res,
                None => {}
            });
    }
    {
        let reads: Vec<_> = ds.encoded_reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
        graph.remove_lightweight_edges(2);
        let cov = ds.coverage.unwrap();
        let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
        graph.remove_zero_copy_elements(cov, &lens, 0.3);
        graph.transitive_edge_reduction();
        let squish = graph.squish_bubbles(15);
        ds.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .for_each(|n| match squish.get(&(n.unit, n.cluster)) {
                Some(res) => n.cluster = *res,
                None => {}
            });
    }
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // println!("unit\tcluster\tto");
    // for ((unit, cluster), to) in squish {
    //     println!("{}\t{}\t{}", unit, cluster, to);
    // }
    Ok(())
}

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
    use haplotyper::assemble::AssembleConfig;
    let config = AssembleConfig::new(2, 100, false, true);
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    graph.remove_lightweight_edges(2);
    let cov = ds.coverage.unwrap();
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    graph.assign_copy_number(cov, &lens);
    // We want to enumerate all the edges in the simple path with copy-numbers more than 1.
    // use haplotyper::Assemble;
    // ds.squish_small_contig(&config);
    // let gfa = ds.assemble(&config);
    // println!("{}", gfa);
    Ok(())
}

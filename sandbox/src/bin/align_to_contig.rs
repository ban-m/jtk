use definitions::*;
use haplotyper::assemble::*;
use haplotyper::model_tune::update_model;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256Plus;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let to_erase: bool = args.get(2).map(|x| x == "true").unwrap_or(false);
    if ds.model_param.is_none() {
        update_model(&mut ds);
    }
    if to_erase {
        ds.encoded_reads
            .iter_mut()
            .flat_map(|n| n.nodes.iter_mut())
            .for_each(|n| n.cluster = 0);
    }
    let config = AssembleConfig::new(1, 100, false, false, 3, 4f64);
    align_to_draft(&ds, &config);
    Ok(())
}
pub fn align_to_draft(ds: &DataSet, c: &AssembleConfig) {
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let mut graph = DitchGraph::new(&reads, &ds.selected_chunks, ds.read_type, c);
    graph.remove_lightweight_edges(2, true);
    let cov = ds.coverage.unwrap_or(30.0);
    let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(4395);
    graph.assign_copy_number_mst(cov, &mut rng);
    // graph.clean_up_graph_for_assemble(cov, &reads, c, ds.read_type);
    eprintln!("{graph}");
    eprintln!("CC:{}", graph.cc());
    assert!(graph.sanity_check());
    let (mut segments, _, _, _, mut encs) = graph.spell(c);
    segments
        .iter()
        .zip(encs.iter())
        .for_each(|(s, e)| assert_eq!(s.sid, e.id));
    while segments.len() > 1 {
        segments.pop();
        encs.pop();
    }
    let config = haplotyper::assemble::consensus::PolishConfig::default();
    haplotyper::assemble::consensus::polish_segment(&ds, &segments, &encs, &config);
}

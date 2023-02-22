use definitions::DataSet;
use haplotyper::model_tune::ModelFit;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // haplotyper::misc::update_coverage(&mut ds);
    // ds.update_models_on_both_strands();
    // println!("{}", serde_json::ser::to_string(&ds).unwrap());
    let hmm = ds.get_model_on_both_strands();
    // let cid = 0;
    // let ref_chunk = ds.selected_chunks.iter().find(|c| c.id == cid).unwrap();
    for ref_chunk in ds.selected_chunks.iter() {
        let cid = ref_chunk.id;
        let (mut reads, mut strands, mut ops) = (vec![], vec![], vec![]);
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            if node.chunk == cid {
                reads.push(node.seq());
                strands.push(node.is_forward);
                ops.push(haplotyper::misc::ops_to_kiley(&node.cigar));
            }
        }
        let config = kiley::hmm::HMMPolishConfig::new(30, reads.len(), 0);
        let start = std::time::Instant::now();
        hmm.polish_until_converge_antidiagonal(
            ref_chunk.seq(),
            &reads,
            &mut ops,
            &strands,
            &config,
        );
        let end = std::time::Instant::now();
        eprintln!("{cid}\t{}", (end - start).as_millis());
    }
    Ok(())
}

use definitions::DataSet;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::repeat_masking::RepeatMask;
    let repetitive_kmer = ds.get_repetitive_kmer();
    for chunk in ds.selected_chunks.iter() {
        let repetitiveness = repetitive_kmer.repetitiveness(chunk.seq());
        println!("{}\t{}", chunk.id, repetitiveness);
    }
    Ok(())
}

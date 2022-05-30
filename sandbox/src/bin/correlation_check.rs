use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::squish_erroneous_clusters::*;
    let config = SquishConfig::default();
    let classify = classify_units(&ds, &config);
    println!("TYPE\tID\tType\tNumCl\tNumCp\tScore");
    for unit in ds.selected_chunks.iter() {
        let id = unit.id;
        let (cl, cp, score) = (unit.cluster_num, unit.copy_num, unit.score);
        let class = classify[&unit.id];
        println!("TYPE\t{id}\t{class}\t{cl}\t{cp}\t{score}");
    }
    Ok(())
}

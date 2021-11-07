use definitions::*;
use std::collections::HashMap;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    let answer: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.contains("hapA")))
        .collect();
    use haplotyper::hapcut::HapCut;
    use haplotyper::hapcut::HapCutConfig;
    let mut config = HapCutConfig::default();
    for seed in vec![342094, 4234, 56, 43290] {
        println!("=========");
        config.set_seed(seed);
        let consis = ds.hapcut(&config);
        println!("SCORE\t{}", consis);
        let mut result: HashMap<_, [u32; 2]> = HashMap::new();
        for &Assignment { id, cluster } in ds.assignments.iter() {
            *result
                .entry(cluster)
                .or_default()
                .get_mut(answer[&id] as usize)
                .unwrap() += 1;
        }
        for (cluster, haps) in result.iter() {
            println!("{}\t{}\t{}", cluster, haps[0], haps[1]);
        }
    }
    Ok(())
}

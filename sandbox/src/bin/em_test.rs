use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::em_correction::*;
    let repeat_num = 1;
    let coverage_thr = 3;
    use log::*;
    debug!("RESULT\tK\tLK\tCov");
    for unit_id in vec![1438, 1494, 24] {
        debug!("{}", unit_id);
        let ref_unit = ds.selected_chunks.iter().find(|n| n.id == unit_id).unwrap();
        let reads: Vec<_> = ds
            .encoded_reads
            .iter()
            .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
            .collect();
        let k = ref_unit.cluster_num;
        let (new_clustering, lk, cluster_num) = (1..=k)
            .flat_map(|k| std::iter::repeat(k).take(repeat_num))
            .map(|k| {
                let seed = unit_id * k as u64;
                let config = Config::new(repeat_num, seed, k, unit_id, coverage_thr);
                em_clustering(&reads, &config)
            })
            .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
            .unwrap();
        let mean_lk = lk / new_clustering.len() as f64;
        debug!(
            "RESULT\t{}\t{}\t{}",
            cluster_num,
            mean_lk,
            new_clustering.len()
        );
    }
    Ok(())
}

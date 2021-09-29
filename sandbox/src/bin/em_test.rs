use definitions::*;
use haplotyper::ClusteringCorrection;
// use rand::SeedableRng;
// use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    let mut pvalues = haplotyper::unit_correlation::calc_p_values(&ds, 5);
    pvalues.retain(|_, pvalue| 0.01 < *pvalue);
    let to_squish = pvalues;
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .filter(|n| to_squish.contains_key(&n.unit))
        .for_each(|n| n.cluster = 0);
    let target: HashSet<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    ds = ds.correct_clustering_em_on_selected(10, 5, true, &target);
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    let mut result: HashMap<_, [usize; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let is_hap1 = id2desc[&read.id].contains("255v2") as usize;
        for node in read.nodes.iter().filter(|n| target.contains(&n.unit)) {
            result.entry((node.unit, node.cluster)).or_default()[is_hap1] += 1;
        }
    }
    println!("unit\tcluster\thap1\thap2");
    for ((unit, cluster), counts) in result {
        println!("{}\t{}\t{}\t{}", unit, cluster, counts[0], counts[1]);
    }
    Ok(())
}

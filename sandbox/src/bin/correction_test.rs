const IS_MOCK: bool = false;
use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        assert!(node.posterior.iter().all(|x| !x.is_nan()));
    }
    println!("UNIT\tid\tunit\tcluster\thap1\thap2\tpurity\tscore");
    {
        let start = std::time::Instant::now();
        let mut ds = ds.clone();
        use haplotyper::dirichlet_mixture::{ClusteringConfig, DirichletMixtureCorrection};
        let config = ClusteringConfig::new(5, 10, 5);
        // use std::collections::HashSet;
        // let selection: HashSet<_> = vec![574].into_iter().collect();
        // ds.correct_clustering_on_selected(&config, &selection);
        ds.correct_clustering(&config);
        let end = std::time::Instant::now();
        eprintln!("DIRTIME\t{}", (end - start).as_secs_f32());
        dump(&ds, "dir_mixture_1");
    }
    // {
    //     let mut ds = ds.clone();
    //     use haplotyper::dirichlet_mixture_old::{ClusteringConfig, DirichletMixtureCorrection};
    //     let config = ClusteringConfig::new(5, 10, 5);
    //     ds.correct_clustering(&config);
    //     dump(&ds, "dir_mixture_2");
    // }
    dump(&ds, "before");
    Ok(())
}

fn dump(ds: &DataSet, id: &str) {
    use std::collections::HashMap;
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
    for ((unit, cluster), counts) in counts.iter() {
        let score = score[unit];
        let total = counts[0] + counts[1];
        let pur = counts[0].max(counts[1]) as f64 / total as f64;
        println!(
            "UNIT\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            id, unit, cluster, counts[0], counts[1], pur, score,
        );
    }
}

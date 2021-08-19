use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::local_clustering;
    let units: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    // let id2desc: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|read| (read.id, &read.desc))
    //     .collect();
    for unit in units {
        let unit = ds.selected_chunks.iter().find(|u| u.id == unit).unwrap();
        let mut config =
            local_clustering::kmeans::ClusteringConfig::new(100, unit.cluster_num as u8, 28f64);
        use rand::SeedableRng;
        use rand_xoshiro::Xoroshiro128Plus;
        let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(293);
        let (queries, ids): (Vec<_>, Vec<_>) = ds
            .encoded_reads
            .iter()
            .flat_map(|r| {
                r.nodes
                    .iter()
                    .filter(|n| n.unit == unit.id)
                    .map(|n| (n.seq(), r.id))
                    .collect::<Vec<_>>()
            })
            .unzip();
        log::debug!("{}", queries.len());
        let (asn, _, _) =
            local_clustering::kmeans::clustering(&queries, &mut rng, &mut config).unwrap();
        // for (asn, id) in asn.iter().zip(ids) {
        //     let is_hapa = id2desc[&id].contains("252v2") as u8;
        //     println!("{}\t{}", asn, is_hapa);
        // }
    }
    Ok(())
}

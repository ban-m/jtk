// use log::*;
// use poa_hmm::*;
use definitions::*;
use log::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    debug!("{:?}", end - start);
    let id2name: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.to_string()))
        .collect();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(0394380);
    for unit in ds.selected_chunks.iter() {
        let (mut seq, mut ids) = (vec![], vec![]);
        for read in ds.encoded_reads.iter() {
            for node in read.nodes.iter().filter(|n| n.unit == unit.id) {
                seq.push(node.seq());
                ids.push(read.id);
            }
        }
        use haplotyper::local_clustering::unit_clustering_ccs_kmervec;
        let initial: Vec<u8> = seq
            .iter()
            .map(|_| rng.gen::<u8>() % unit.cluster_num as u8)
            .collect();
        let asn = unit_clustering_ccs_kmervec(&seq, unit.cluster_num as u8, 8, &initial);
        for (i, (id, asn)) in ids.iter().zip(asn.iter()).enumerate() {
            println!("{}\t{}\t{}\t{}", unit.id, i, id2name[id], asn);
        }
    }
    Ok(())
}

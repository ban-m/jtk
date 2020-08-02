use haplotyper::clustering_by_kmeans;
use haplotyper::eread::*;
use haplotyper::ClusteringConfig;
use poa_hmm::*;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(322321);
    let chain_len = 20;
    let ref_unit = 0;
    let mut c = ClusteringConfig::default();
    // c.read_type = haplotyper::ReadType::CLR;
    c.read_type = haplotyper::ReadType::CCS;
    c.cluster_num = 2;
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, 100))
        .collect();
    let template2: Vec<_> = template
        .iter()
        .enumerate()
        .map(|(idx, x)| {
            if idx == 5 {
                gen_sample::introduce_errors(x, &mut rng, 1, 0, 0)
            } else {
                x.clone()
            }
        })
        .collect();
    let template3 = template.clone();
    // let profile = &gen_sample::PROFILE;
    let profile = &gen_sample::CCS_PROFILE;
    let seq1: Vec<_> = (0..30)
        .map(|_| {
            let chunks: Vec<_> = template
                .iter()
                .enumerate()
                .map(|(pos, seq)| {
                    let seq = gen_sample::introduce_randomness(seq, &mut rng, profile);
                    Chunk { pos, seq }
                })
                .collect();
            let cluster = 0;
            ChunkedUnit { cluster, chunks }
        })
        .collect();
    let seq2: Vec<_> = (0..30)
        .map(|_| {
            let chunks: Vec<_> = template2
                .iter()
                .enumerate()
                .map(|(pos, seq)| {
                    let seq = gen_sample::introduce_randomness(seq, &mut rng, profile);
                    Chunk { pos, seq }
                })
                .collect();
            let cluster = 0;
            ChunkedUnit { cluster, chunks }
        })
        .collect();
    let seq3: Vec<_> = (0..30)
        .map(|_| {
            let chunks: Vec<_> = template3
                .iter()
                .enumerate()
                .map(|(pos, seq)| {
                    let seq = gen_sample::introduce_randomness(seq, &mut rng, profile);
                    Chunk { pos, seq }
                })
                .collect();
            let cluster = 0;
            ChunkedUnit { cluster, chunks }
        })
        .collect();
    c.seed = 121;
    let mut data = vec![seq1.clone(), seq2.clone()].concat();
    clustering_by_kmeans(&mut data, chain_len, &c, ref_unit);
    for (idx, d) in data.iter().enumerate() {
        eprintln!("{}\t{}", idx, d.cluster);
    }
    let mut data = vec![seq1.clone(), seq3.clone()].concat();
    clustering_by_kmeans(&mut data, chain_len, &c, ref_unit);
    for (idx, d) in data.iter().enumerate() {
        eprintln!("{}\t{}", idx, d.cluster);
    }
    Ok(())
}

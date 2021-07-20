use log::debug;
use poa_hmm::gen_sample;
use rand::Rng;
pub fn generate_mul_data<T: Rng>(
    templates: &[Vec<Vec<u8>>],
    test_num: usize,
    rng: &mut T,
    _probs: &[f64],
    profile: &poa_hmm::gen_sample::Profile,
) -> (Vec<haplotyper::eread::ChunkedUnit>, Vec<u8>) {
    let choices: Vec<_> = (0..templates.len()).collect();
    use rand::seq::SliceRandom;
    let mut answer: Vec<_> = (0..test_num)
        .filter_map(|_| choices.choose(rng))
        .map(|&x| x as u8)
        .collect();
    answer.sort();
    let mut gen = |t: &[Vec<u8>]| {
        t.iter()
            .map(|e| gen_sample::introduce_randomness(e, rng, profile))
            .enumerate()
            .map(|(pos, seq)| haplotyper::eread::Chunk { pos, seq })
            .collect::<Vec<_>>()
    };
    let dataset: Vec<_> = answer
        .iter()
        .map(|&idx| gen(&templates[idx as usize]))
        .map(|chunks| haplotyper::eread::ChunkedUnit { chunks, cluster: 0 })
        .collect();
    assert_eq!(dataset.len(), answer.len());
    debug!("Index1\tIndex2\tDist");
    let clusters = templates.len();
    (0..clusters).for_each(|i| {
        (i + 1..clusters).for_each(|j| {
            let dist = templates[i]
                .iter()
                .zip(templates[j].iter())
                .enumerate()
                .map(|(idx, (t1, t2))| {
                    let d = bio_utils::alignments::edit_dist(t1, t2);
                    if d != 0 {
                        debug!("VARPOS\t{}\t{}", idx, d);
                        // debug!("{}", String::from_utf8_lossy(t1));
                        // debug!("{}", String::from_utf8_lossy(t2));
                    }
                    d
                })
                .sum::<u32>();
            debug!("{}\t{}\t{}", i, j, dist);
        });
    });
    (dataset, answer)
}

pub fn clustering(reads: &[(usize, Vec<Vec<u8>>)], chain: usize, k: usize, s: u64) -> Vec<usize> {
    let mut data: Vec<_> = reads
        .iter()
        .map(|r| {
            use haplotyper::local_clustering::eread::*;
            let chunks: Vec<_> =
                r.1.iter()
                    .enumerate()
                    .map(|(i, c)| Chunk {
                        pos: i,
                        seq: c.to_vec(),
                    })
                    .collect();
            ChunkedUnit { cluster: 0, chunks }
        })
        .collect();
    let unit = definitions::Unit {
        id: 0,
        seq: String::new(),
        cluster_num: k,
    };
    // haplotyper::local_clustering::initial_clustering(&mut data, &unit, (k, chain), s);
    use haplotyper::local_clustering::clustering_by_kmeans_em;
    let mut c = haplotyper::local_clustering::ClusteringConfig::default();
    c.poa_config = poa_hmm::DEFAULT_CONFIG;
    c.read_type = definitions::ReadType::CLR;
    c.variant_num = 2;
    c.cluster_num = k;
    clustering_by_kmeans_em(&mut data, chain, &c, &unit, s);
    data.iter().map(|d| d.cluster).collect()
}

pub fn em_clustering(
    xs: &[Vec<usize>],
    trials: usize,
    k: usize,
    cluster: usize,
    s: u64,
) -> (Vec<usize>, f64) {
    haplotyper::local_clustering::em_clustering(xs, cluster, (s, trials, k))
}

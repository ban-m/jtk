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
                        debug!("{}\t{}", idx, d);
                    }
                    d
                })
                .sum::<u32>();
            debug!("{}\t{}\t{}", i, j, dist);
        });
    });
    (dataset, answer)
}

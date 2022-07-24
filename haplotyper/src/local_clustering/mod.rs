use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
mod config;
pub use config::*;
pub mod kmeans;
pub mod normalize;
/// Return rand index.
pub fn rand_index(label: &[usize], pred: &[usize]) -> f64 {
    assert_eq!(label.len(), pred.len());
    let mut both_same_pair = 0;
    let mut both_diff_pair = 0;
    for (i, (l1, p1)) in label.iter().zip(pred.iter()).enumerate() {
        for (l2, p2) in label.iter().zip(pred.iter()).take(i) {
            if l1 == l2 && p1 == p2 {
                both_same_pair += 1;
            } else if l1 != l2 && p1 != p2 {
                both_diff_pair += 1;
            }
        }
    }
    let len = label.len();
    (both_same_pair + both_diff_pair) as f64 / (len * (len - 1) / 2) as f64
}

pub fn adjusted_rand_index(label: &[usize], pred: &[usize]) -> f64 {
    assert_eq!(label.len(), pred.len());
    let lab_max = *label.iter().max().unwrap();
    let pred_max = *pred.iter().max().unwrap();
    let mut cont_table = vec![vec![0; pred_max + 1]; lab_max + 1];
    let mut lab_sum = vec![0; lab_max + 1];
    let mut pred_sum = vec![0; pred_max + 1];

    for (&lab, &pred) in label.iter().zip(pred.iter()) {
        cont_table[lab][pred] += 1;
        lab_sum[lab] += 1;
        pred_sum[pred] += 1;
    }
    fn choose(x: &usize) -> usize {
        (x.max(&1) - 1) * x / 2
    }
    let lab_match: usize = lab_sum.iter().map(choose).sum();
    let pred_match: usize = pred_sum.iter().map(choose).sum();
    let num_of_pairs = choose(&label.len());
    let both_match: usize = cont_table.iter().flatten().map(choose).sum();
    let denom = num_of_pairs * (lab_match + pred_match) / 2 - lab_match * pred_match;
    let numer = num_of_pairs * both_match - lab_match * pred_match;
    numer as f64 / denom as f64
}

pub trait LocalClustering {
    fn local_clustering(&mut self);
}

impl LocalClustering for DataSet {
    fn local_clustering(&mut self) {
        let selection: HashSet<_> = self.selected_chunks.iter().map(|x| x.id).collect();
        local_clustering_selected(self, &selection);
    }
}

fn set_coverage(ds: &mut DataSet) {
    let mut counts: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    let cov = {
        let mut counts: Vec<_> = counts.values().copied().collect();
        counts.sort_unstable();
        counts[counts.len() / 2] as f64 / 2f64
    };
    debug!("LOCALCLUSTERING\tSetCoverage\t{cov}");
    ds.coverage = Some(cov);
}

/// Selection: HashSet of the chunk ID to be clustered on.
pub fn local_clustering_selected(ds: &mut DataSet, selection: &HashSet<u64>) {
    if selection.is_empty() {
        return;
    }
    if ds.coverage.is_none() {
        set_coverage(ds);
    }
    if ds.model_param.is_none() {
        crate::model_tune::update_model(ds);
    }
    let hmm = crate::model_tune::get_model(ds).unwrap();
    // use crate::stats::Stats;
    let gain = estimate_minimum_gain(&hmm); //ds.error_rate()
    debug!("MinGain\t{:.3}", gain);
    let mut pileups: HashMap<u64, Vec<&mut Node>> =
        selection.iter().map(|&id| (id, vec![])).collect();
    let chunks: HashMap<u64, _> = ds
        .selected_chunks
        .iter()
        .filter(|c| selection.contains(&c.id))
        .map(|c| (c.id, c))
        .collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some(bucket) = pileups.get_mut(&node.unit) {
            bucket.push(node);
        }
    }
    pileups.iter_mut().for_each(|(unit_id, nodes)| {
        nodes.sort_by_cached_key(|node| {
            let ref_unit = chunks.get(unit_id).unwrap();
            let (_, aln, _) = node.recover(ref_unit);
            aln.iter().filter(|&&x| x != b'|').count()
        });
    });
    let coverage = ds.coverage.unwrap();
    let read_type = ds.read_type;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .filter(|(_, units)| !units.is_empty())
        .map(|(&unit_id, units)| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 3490);
            let (seqs, mut ops): (Vec<_>, Vec<_>) = units
                .iter()
                .map(|node| (node.seq(), crate::misc::ops_to_kiley(&node.cigar)))
                .unzip();
            let band_width = read_type.band_width(ref_unit.seq().len());
            let start = std::time::Instant::now();
            use kmeans::ClusteringConfig;
            let copy_num = ref_unit.copy_num as u8;
            let refseq = ref_unit.seq();
            use kiley::bialignment::guided::polish_until_converge_with;
            let cons = polish_until_converge_with(refseq, &seqs, &mut ops, band_width);
            let cons = hmm.polish_until_converge_with(&cons, &seqs, &mut ops, band_width);
            let config = ClusteringConfig::new(band_width / 2, copy_num, coverage, gain);
            use kmeans::*;
            let polished = std::time::Instant::now();
            let strands: Vec<_> = units.iter().map(|n| n.is_forward).collect();
            let (asn, pss, score, k) = if 1 < ref_unit.copy_num {
                let cls = clustering_dev(&cons, &seqs, &mut ops, &strands, &mut rng, &hmm, &config);
                cls.unwrap_or_else(|| panic!("RECORD\t{}\tMISS", unit_id))
            } else {
                (vec![0; units.len()], vec![vec![0f64]; units.len()], 0f64, 1)
            };
            for (node, ps) in units.iter_mut().zip(pss) {
                assert_eq!(ps.len(), k);
                node.posterior = ps;
            }
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            for (node, ops) in units.iter_mut().zip(ops) {
                node.cigar = crate::encode::compress_kiley_ops(&ops).into();
            }
            let end = std::time::Instant::now();
            let polished_time = (polished - start).as_millis();
            let elapsed = (end - start).as_millis();
            let len = cons.len();
            let cov = units.len();
            debug!("RECORD\t{unit_id}\t{elapsed}\t{polished_time}\t{len}\t{score:.3}\t{cov}",);
            (unit_id, (cons, score, k))
        })
        .collect();
    debug!("LC\t{}", consensus_and_clusternum.len());
    for unit in ds.selected_chunks.iter_mut() {
        if let Some(&(ref consensus, score, cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = consensus.to_vec().into();
            unit.score = score;
            unit.cluster_num = cluster_num;
        }
    }
    normalize::normalize_local_clustering(ds);
}

pub fn estimate_minimum_gain(hmm: &kiley::hmm::guided::PairHiddenMarkovModel) -> f64 {
    const SEED: u64 = 23908;
    const SAMPLE_NUM: usize = 1000;
    const SEQ_NUM: usize = 500;
    const LEN: usize = 100;
    const BAND: usize = 25;
    // const FRAC: f64 = 0.01;
    const MIN_REQ: f64 = 1f64;
    let mut medians: Vec<_> = (0..SAMPLE_NUM)
        .into_par_iter()
        .map(|seed| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(SEED + seed as u64);
            let hap1 = kiley::gen_seq::generate_seq(&mut rng, LEN);
            let hap2 = match seed % 2 == 0 {
                true => kiley::gen_seq::introduce_errors(&hap1, &mut rng, 0, 0, 1),
                false => kiley::gen_seq::introduce_errors(&hap1, &mut rng, 0, 1, 0),
            };
            use kiley::gen_seq::Generate;
            let mut lks: Vec<_> = (0..SEQ_NUM)
                .map(|_| {
                    let read = hmm.gen(&hap1, &mut rng);
                    let lk_base = hmm.likelihood(&hap1, &read, BAND);
                    let lk_diff = hmm.likelihood(&hap2, &read, BAND);
                    lk_base - lk_diff
                })
                .collect();
            *lks.select_nth_unstable_by(SEQ_NUM / 2, |x, y| x.partial_cmp(y).unwrap())
                .1
        })
        .collect();
    medians.sort_by(|x, y| x.partial_cmp(y).unwrap());
    debug!("MIN_GAIN\t{:?}", &medians[..6]);
    (medians[2]).max(MIN_REQ)
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn rand_index_test() {
        let pred = [0, 0, 0, 1, 1, 1];
        let answ = [0, 0, 1, 1, 2, 2];
        assert!((0.6666 - rand_index(&pred, &answ)).abs() < 0.0001);
    }
}

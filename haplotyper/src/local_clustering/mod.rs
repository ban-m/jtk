use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
mod config;
pub use config::*;
pub mod kmeans;

/// Return rand index.
pub fn rand_index(label: &[u8], pred: &[u8]) -> f64 {
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

pub trait LocalClustering {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        self,
        c: &ClusteringConfig<F>,
    ) -> Self;
}

impl LocalClustering for DataSet {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        mut self,
        _c: &ClusteringConfig<F>,
    ) -> Self {
        let selection: HashSet<_> = self.selected_chunks.iter().map(|x| x.id).collect();
        local_clustering_selected(&mut self, &selection);
        self
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
    debug!("COVERAGE\t{}\tHAPLOID", cov);
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
    let coverage = ds.coverage.unwrap();
    // Maybe we can train pair-HMM here.
    let hmm = {
        let mut covs: Vec<_> = pileups.iter().map(|x| x.1.len()).collect();
        let (_, cov, _) = covs.select_nth_unstable(pileups.len() / 2);
        let (unit_id, units) = pileups.iter().find(|(_, us)| us.len() == *cov).unwrap();
        debug!("LOCAL\tSAMPLE\t{}\t{}", unit_id, units.len());
        let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
        let ref_unit = chunks.get(&unit_id).unwrap();
        let mut hmm = kiley::gphmm::GPHMM::<Cond>::clr();
        let band_width = 200;
        use kiley::gphmm::*;
        for _ in 0..2 {
            let consensus = take_consensus(ref_unit, &seqs, &hmm);
            hmm = hmm.fit_banded(&consensus, &seqs, band_width);
        }
        hmm
    };
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .map(|(&unit_id, units)| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 25);
            let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
            let cov = seqs.len();
            let ref_unit = chunks.get(&unit_id).unwrap();
            let start = std::time::Instant::now();
            let config = kmeans::ClusteringConfig::new(100, ref_unit.cluster_num as u8, coverage);
            // Maybe it is better to use the original alignment, right?
            let consensus = take_consensus(ref_unit, &seqs, &hmm);
            let (asn, score) = if 1 < ref_unit.cluster_num {
                kmeans::clustering_with_template(&consensus, &seqs, &mut rng, &hmm, &config)
                    .unwrap_or_else(|| panic!("RECORD\t{}", unit_id))
            } else {
                (vec![0; units.len()], 0f64)
            };
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            let end = std::time::Instant::now();
            let elapsed = (end - start).as_secs();
            let len = consensus.len();
            debug!(
                "RECORD\t{}\t{}\t{}\t{:.3}\t{}",
                unit_id, elapsed, len, score, cov
            );
            (unit_id, (consensus, score, config.cluster_num))
        })
        .collect();
    for unit in ds.selected_chunks.iter_mut() {
        if let Some((consensus, score, _cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = String::from_utf8(consensus.to_vec()).unwrap();
            unit.score = *score;
        }
    }
    // NOTE that we do not need to remove the very different reads, as
    // it is already "different" with respect to the clustering information.
    re_encode_reads(ds, &consensus_and_clusternum);
}

pub fn take_consensus<T: std::borrow::Borrow<[u8]>>(
    unit: &Unit,
    reads: &[T],
    hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
) -> Vec<u8> {
    use kiley::polish_chunk_by_parts;
    use kiley::PolishConfig;
    let config = PolishConfig::with_model(100, 0, reads.len(), unit.id, 0, hmm.clone());
    let max_len = reads
        .iter()
        .map(|x| x.borrow().len())
        .max()
        .unwrap_or_else(|| panic!("{},{}", line!(), unit.id));
    match 200 < max_len {
        true => polish_chunk_by_parts(unit.seq(), reads, &config),
        false => unit.seq().to_vec(),
    }
}

fn re_encode_reads(ds: &mut DataSet, consensus: &HashMap<u64, (Vec<u8>, f64, u8)>) {
    ds.encoded_reads
        .par_iter_mut()
        .flat_map(|r| r.nodes.par_iter_mut())
        .for_each(|node| {
            let cons = match consensus.get(&node.unit) {
                Some((cons, _, _)) => cons,
                None => return,
            };
            let band_size = (cons.len() / 10).max(5);
            let ops =
                kiley::bialignment::global_banded(cons, node.seq(), 2, -5, -6, -1, band_size).1;
            // I think we need not to filtering out weak alignment node,
            // it is task to the downstream algorithm.
            node.cigar = crate::encode::compress_kiley_ops(&ops);
        });
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

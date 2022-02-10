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
    let hmm = estimate_model_parameters(ds.read_type, &pileups, &chunks);
    let band_width = ds.read_type.band_width();
    let read_type = ds.read_type;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .filter(|(_, units)| !units.is_empty())
        .map(|(&unit_id, units)| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 25);
            let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
            let cov = seqs.len();
            let ref_unit = chunks.get(&unit_id).unwrap();
            let start = std::time::Instant::now();
            use kmeans::ClusteringConfig;
            let copy_num = ref_unit.copy_num as u8;
            let config = ClusteringConfig::new(band_width, copy_num, coverage, read_type);
            let consensus = {
                let mut seqs_with_diff: Vec<_> = units
                    .iter()
                    .map(|node| {
                        let (_, aln, _) = node.recover(ref_unit);
                        let diff: usize = aln.iter().filter(|&&x| x != b'|').count();
                        (node.seq(), diff)
                    })
                    .collect();
                seqs_with_diff.sort_by_key(|x| x.1);
                let seqs: Vec<_> = seqs_with_diff.iter().take(50).map(|x| x.0).collect();
                take_consensus(ref_unit, &seqs, band_width, &hmm)
            };
            let (asn, pss, score, k) = if 1 < ref_unit.copy_num {
                kmeans::clustering_dev(&consensus, &seqs, &mut rng, &hmm, &config)
                    .unwrap_or_else(|| panic!("RECORD\t{}", unit_id))
            } else {
                (vec![0; units.len()], vec![vec![0f64]; units.len()], 0f64, 1)
            };
            for (node, ps) in units.iter_mut().zip(pss) {
                node.posterior = ps;
            }
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
            (unit_id, (consensus, score, k))
        })
        .collect();
    for unit in ds.selected_chunks.iter_mut() {
        if let Some((consensus, score, cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = String::from_utf8(consensus.to_vec()).unwrap();
            unit.score = *score;
            unit.cluster_num = *cluster_num as usize;
        }
    }
    // NOTE that we do not need to remove the very different reads, as
    // it is already "different" with respect to the clustering information.
    re_encode_reads(ds, &consensus_and_clusternum);
    // Normalize local clustering.
    normalize::normalize_local_clustering(ds);
}

use kiley::gphmm::Cond;
pub fn estimate_model_parameters<N: std::borrow::Borrow<Node>>(
    read_type: ReadType,
    pileups: &HashMap<u64, Vec<N>>,
    chunks: &HashMap<u64, &Unit>,
) -> kiley::gphmm::GPHMM<Cond> {
    let mut covs: Vec<_> = pileups.iter().map(|x| x.1.len()).collect();
    let (_, &mut cov, _) = covs.select_nth_unstable(pileups.len() / 2);
    let seqs_and_ref_units: Vec<_> = pileups
        .iter()
        .filter(|(_, us)| (cov.max(1) - 1..cov + 2).contains(&us.len()))
        .map(|(id, us)| {
            let seqs: Vec<_> = us.iter().map(|n| n.borrow().seq()).take(100).collect();
            let ref_unit = chunks.get(id).unwrap();
            (ref_unit, seqs)
        })
        .take(2)
        .collect();
    for (chunk, units) in seqs_and_ref_units.iter() {
        debug!("LOCAL\tSAMPLE\t{}\t{}", chunk.id, units.len());
    }
    use kiley::gphmm::*;
    let mut hmm = match read_type {
        ReadType::CCS => kiley::gphmm::GPHMM::<Cond>::ccs(),
        ReadType::None | ReadType::CLR => kiley::gphmm::GPHMM::<Cond>::clr(),
        ReadType::ONT => kiley::gphmm::GPHMM::<Cond>::ont(),
    };
    for _ in 0..2 {
        for (ref_unit, seqs) in seqs_and_ref_units.iter() {
            let band_width = read_type.band_width() * 2;
            let consensus = take_consensus(ref_unit, seqs, band_width, &hmm);
            hmm = hmm.fit_banded(&consensus, seqs, band_width);
        }
    }
    debug!("HMM\t{}", hmm);
    hmm
}

pub fn take_consensus<T: std::borrow::Borrow<[u8]>>(
    unit: &Unit,
    reads: &[T],
    band_width: usize,
    hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
) -> Vec<u8> {
    use kiley::polish_chunk_by_parts;
    use kiley::PolishConfig;
    let config = PolishConfig::with_model(band_width, 0, reads.len(), unit.id, 0, hmm.clone());
    let max_len = reads
        .iter()
        .map(|x| x.borrow().len())
        .max()
        .unwrap_or_else(|| panic!("{},{}", unit.seq, unit.id));
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

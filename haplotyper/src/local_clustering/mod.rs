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
    pileups.iter_mut().for_each(|(unit_id, nodes)| {
        nodes.sort_by_cached_key(|node| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            let (_, aln, _) = node.recover(ref_unit);
            aln.iter().filter(|&&x| x != b'|').count()
        });
    });
    let coverage = ds.coverage.unwrap();
    let hmm = estimate_model_parameters_neo(ds.read_type, &pileups, &chunks);
    let read_type = ds.read_type;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .filter(|(_, units)| !units.is_empty())
        .map(|(&unit_id, units)| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 25);
            let (seqs, mut ops): (Vec<_>, Vec<_>) = units
                .iter()
                .map(|node| {
                    let mut ops = vec![];
                    for op in node.cigar.iter() {
                        match op {
                            Op::Match(l) => {
                                ops.extend(std::iter::repeat(kiley::Op::Match).take(*l))
                            }
                            Op::Del(l) => ops.extend(std::iter::repeat(kiley::Op::Del).take(*l)),
                            Op::Ins(l) => ops.extend(std::iter::repeat(kiley::Op::Ins).take(*l)),
                        }
                    }
                    (node.seq(), ops)
                })
                .unzip();
            let cov = seqs.len();
            let band_width = read_type.band_width(ref_unit.seq().len()) / 2;
            let start = std::time::Instant::now();
            use kmeans::ClusteringConfig;
            let copy_num = ref_unit.copy_num as u8;
            let consensus =
                hmm.polish_until_converge_with(ref_unit.seq(), &seqs, &mut ops, band_width);
            let config = ClusteringConfig::new(band_width / 2, copy_num, coverage, read_type);
            let (asn, pss, score, k) = if 1 < ref_unit.copy_num {
                kmeans::clustering_neo(&consensus, &seqs, &mut ops, &mut rng, &hmm, &config)
                    .unwrap_or_else(|| panic!("RECORD\t{}\tMISS", unit_id))
            } else {
                (vec![0; units.len()], vec![vec![0f64]; units.len()], 0f64, 1)
            };
            for (node, ps) in units.iter_mut().zip(pss) {
                node.posterior = ps;
            }
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            // Change cigar string...
            let end = std::time::Instant::now();
            let elapsed = (end - start).as_millis();
            let len = consensus.len();
            debug!(
                "RECORD\t{}\t{}\t{}\t{:.3}\t{}",
                unit_id, elapsed, len, score, cov
            );
            (unit_id, (consensus, score, k))
        })
        .collect();
    debug!("LC\t{}", consensus_and_clusternum.len());
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

pub fn estimate_model_parameters_neo<N: std::borrow::Borrow<Node>>(
    read_type: ReadType,
    pileups: &HashMap<u64, Vec<N>>,
    chunks: &HashMap<u64, &Unit>,
) -> kiley::hmm::guided::PairHiddenMarkovModel {
    let mut covs: Vec<_> = pileups.iter().map(|x| x.1.len()).collect();
    let (_, &mut cov, _) = covs.select_nth_unstable(pileups.len() / 2);
    let mut seqs_and_ref_units: Vec<_> = pileups
        .iter()
        .filter(|(_, us)| (cov.max(1) - 1..cov + 2).contains(&us.len()))
        .map(|(id, us)| (chunks.get(id).unwrap(), us))
        .collect();
    seqs_and_ref_units.sort_by_cached_key(|c| c.0.id);
    seqs_and_ref_units.truncate(2);
    for (chunk, units) in seqs_and_ref_units.iter() {
        debug!("LOCAL\tSAMPLE\t{}\t{}", chunk.id, units.len());
    }
    let mut hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let mut polishing_pairs: Vec<_> = seqs_and_ref_units
        .iter()
        .map(|(ref_unit, nodes)| {
            let band_width = 2 * read_type.band_width(ref_unit.seq().len());
            let ops: Vec<Vec<_>> = nodes
                .iter()
                .map(|n| {
                    n.borrow()
                        .cigar
                        .iter()
                        .flat_map(|op| match *op {
                            definitions::Op::Del(l) => std::iter::repeat(kiley::Op::Del).take(l),
                            definitions::Op::Ins(l) => std::iter::repeat(kiley::Op::Ins).take(l),
                            definitions::Op::Match(l) => {
                                std::iter::repeat(kiley::Op::Match).take(l)
                            }
                        })
                        .collect()
                })
                .collect();
            let seqs: Vec<_> = nodes.iter().map(|n| n.borrow().seq()).collect();
            (ref_unit.seq().to_vec(), seqs, ops, band_width)
        })
        .collect();
    polishing_pairs
        .par_iter_mut()
        .for_each(|(consensus, seqs, ops, bw)| {
            *consensus = hmm.polish_until_converge_with(consensus, seqs, ops, *bw);
        });
    debug!("POLISHED");
    for _ in 0..2 {
        for (consensus, seqs, ops, bw) in polishing_pairs.iter_mut() {
            hmm.fit_naive_with_par(consensus, seqs, ops, *bw);
        }
    }
    debug!("HMM\n{}", hmm);
    hmm
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
            let band_width = read_type.band_width(ref_unit.seq().len()) * 2;
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
    let max_len = reads
        .iter()
        .map(|x| x.borrow().len())
        .max()
        .unwrap_or_else(|| panic!("{},{}", unit.seq, unit.id));
    match 200 < max_len {
        true => {
            use kiley::polish_chunk_by_parts;
            use kiley::PolishConfig;
            let config =
                PolishConfig::with_model(band_width, 0, reads.len(), unit.id, 0, hmm.clone());
            polish_chunk_by_parts(unit.seq(), reads, &config)
        }
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

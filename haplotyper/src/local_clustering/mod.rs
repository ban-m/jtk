use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
mod config;
pub use config::*;
pub mod kmeans;
pub mod normalize;

pub trait LocalClustering {
    fn local_clustering(&mut self);
}

impl LocalClustering for DataSet {
    fn local_clustering(&mut self) {
        let selection: HashSet<_> = self.selected_chunks.iter().map(|x| x.id).collect();
        local_clustering_selected(self, &selection);
    }
}

fn pileup_nodes<'a>(
    ds: &'a mut DataSet,
    selection: &HashSet<u64>,
) -> (HashMap<u64, Vec<&'a mut Node>>, HashMap<u64, &'a Unit>) {
    let chunks: HashMap<u64, _> = ds
        .selected_chunks
        .iter()
        .filter(|c| selection.contains(&c.id))
        .map(|c| (c.id, c))
        .collect();
    let mut pileups: HashMap<u64, Vec<&mut Node>> =
        selection.iter().map(|&id| (id, vec![])).collect();

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
    (pileups, chunks)
}

const SEED: u64 = 309423;
const SEQ_LEN: usize = 100;
const BAND: usize = 10;
const HOMOP_LEN: usize = 3;
/// Selection: HashSet of the chunk ID to be clustered on.
pub fn local_clustering_selected(ds: &mut DataSet, selection: &HashSet<u64>) {
    use kmeans::*;
    if selection.is_empty() {
        return;
    }
    crate::misc::update_coverage(ds);
    if ds.model_param.is_none() {
        crate::model_tune::update_model(ds);
    }
    let coverage = ds.coverage.unwrap();
    let read_type = ds.read_type;
    let hmm = crate::model_tune::get_model(ds).unwrap();
    let (mut pileups, chunks) = pileup_nodes(ds, selection);
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
            let copy_num = ref_unit.copy_num as u8;
            let refseq = ref_unit.seq();
            let (cons, hmm) = prep_consensus(&hmm, refseq, &seqs, &mut ops, band_width);
            let gains =
                crate::likelihood_gains::estimate_gain(&hmm, SEED, SEQ_LEN, BAND, HOMOP_LEN);
            let config = ClusteringConfig::new(band_width / 2, copy_num, coverage, &gains);
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
                node.cigar = crate::misc::kiley_op_to_ops(&ops);
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

// TODO: this function is, very very slow. Please fasten this function, please.
type Phmm = kiley::hmm::guided::PairHiddenMarkovModel;
fn prep_consensus(
    hmm: &Phmm,
    draft: &[u8],
    seqs: &[&[u8]],
    ops: &mut [Vec<kiley::Op>],
    band_width: usize,
) -> (Vec<u8>, Phmm) {
    let mut hmm = hmm.clone();
    let mut cons = draft.to_vec();
    for _t in 0..2 {
        hmm.fit_naive_with(&cons, seqs, ops, band_width / 2);
        cons = hmm.polish_until_converge_with(&cons, seqs, ops, band_width / 2);
    }
    (cons, hmm)
}

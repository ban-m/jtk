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
    let (pileups, chunks) = pileup_nodes(ds, selection);
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .into_par_iter()
        .filter(|(_, units)| !units.is_empty())
        .map(|(unit_id, mut units)| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            // let (mut units, mut unused_units) =
            //     partition_by_lk_capability(ref_unit, units, read_type, &hmm);
            assert!(!units.is_empty());
            let (uid, (template, score, cl_num)) =
                clustering_on_pileup(&mut units, ref_unit, read_type, &hmm, coverage);
            // update_by_removed_nodes(&mut unused_units, cl_num, &template);
            (uid, (template, score, cl_num))
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

// fn partition_by_lk_capability<'a>(
//     ref_unit: &Unit,
//     units: Vec<&'a mut Node>,
//     read_type: ReadType,
//     hmm: &Phmm,
// ) -> (Vec<&'a mut Node>, Vec<&'a mut Node>) {
//     let refseq = ref_unit.seq();
//     let band_width = read_type.band_width(ref_unit.seq().len());
//     units.into_iter().partition(|node| {
//         let ops = crate::misc::ops_to_kiley(&node.cigar);
//         let forward = hmm.likelihood_guided(refseq, node.seq(), &ops, band_width / 2);
//         let forward_ng = forward.is_nan() || forward.is_infinite();
//         let backward = hmm.likelihood_guided_post(refseq, node.seq(), &ops, band_width / 2);
//         let backward_ng = backward.is_nan() || backward.is_infinite();
//         !forward_ng && !backward_ng
//     })
// }

const UPPER_COPY_NUM: usize = 8;
fn clustering_on_pileup(
    units: &mut [&mut Node],
    ref_unit: &Unit,
    read_type: ReadType,
    hmm: &Phmm,
    coverage: f64,
) -> (u64, (Vec<u8>, f64, usize)) {
    use kmeans::*;
    let refseq = ref_unit.seq();
    let band_width = read_type.band_width(ref_unit.seq().len());
    let unit_id = ref_unit.id;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 3490);
    let (seqs, mut ops): (Vec<_>, Vec<_>) = units
        .iter()
        .map(|node| (node.seq(), crate::misc::ops_to_kiley(&node.cigar)))
        .unzip();
    let start = std::time::Instant::now();
    let copy_num = ref_unit.copy_num;
    let (cons, hmm) = prep_consensus(&hmm, refseq, &seqs, &mut ops, band_width);
    let polished = std::time::Instant::now();
    let strands: Vec<_> = units.iter().map(|n| n.is_forward).collect();
    let gains = crate::likelihood_gains::estimate_gain(&hmm, SEED, SEQ_LEN, BAND, HOMOP_LEN);
    let per_cluster_cov = match copy_num {
        0 | 1 | 2 => seqs.len() as f64 / copy_num as f64,
        _ => (seqs.len() as f64 / copy_num as f64).max(coverage),
    };
    let config = ClusteringConfig::new(band_width / 2, copy_num, coverage, per_cluster_cov, &gains);
    let (asn, pss, score, k) =
        clustering_recursive(&cons, &seqs, &ops, &strands, &mut rng, &hmm, &config);
    update_by_clusterings(units, &asn, &ops, &pss);
    let end = std::time::Instant::now();
    let polished_time = (polished - start).as_millis();
    let elapsed = (end - start).as_millis();
    let (len, cov) = (cons.len(), units.len());
    debug!("RECORD\t{unit_id}\t{elapsed}\t{polished_time}\t{len}\t{score:.3}\t{cov}",);
    (unit_id, (cons, score, k))
}

type ClusteringDevResult = (Vec<usize>, Vec<Vec<f64>>, f64, usize);
fn clustering_recursive<R: rand::Rng>(
    cons: &[u8],
    seqs: &[&[u8]],
    ops: &[Vec<kiley::Op>],
    strands: &[bool],
    rng: &mut R,
    hmm: &Phmm,
    config: &kmeans::ClusteringConfig,
) -> ClusteringDevResult {
    use kmeans::*;
    if config.copy_num < UPPER_COPY_NUM {
        clustering_dev(cons, seqs, ops, &strands, rng, hmm, config)
    } else {
        trace!("RECURSE");
        const BRANCH_NUM: usize = 4;
        let mut rec_config = config.clone();
        rec_config.copy_num = BRANCH_NUM;
        let band_width = config.band_width;
        let (asn, pss, score, k) = clustering_dev(cons, seqs, ops, &strands, rng, hmm, &rec_config);
        let copy_numbers = estim_copy_num(&asn, k, config.copy_num, config.coverage);
        trace!("RECURSE\tCOPYNUM\t{:?}", copy_numbers);
        if k <= 1 {
            return (asn, pss, score, k);
        }
        let recurred_clusterings: Vec<_> = copy_numbers
            .iter()
            .enumerate()
            .map(|(k, &cp)| {
                let (seqs, mut ops, strands) = filter_sub_clusters(seqs, ops, strands, &asn, k);
                let (cons, hmm) = prep_consensus(hmm, &cons, &seqs, &mut ops, band_width);
                let mut config = config.clone();
                config.copy_num = cp;
                clustering_recursive(&cons, &seqs, &ops, &strands, rng, &hmm, &config)
            })
            .collect();
        let sub_scores: f64 = recurred_clusterings.iter().map(|x| x.2).sum();
        let total_score = sub_scores + score;
        let cluster_nums: Vec<_> = recurred_clusterings.iter().map(|x| x.3).collect();
        trace!("RECURSE\tCLNUM\t{:?}\t{}", cluster_nums, config.copy_num);
        let total_cluster_num: usize = cluster_nums.iter().sum();
        let offsets = cumsum(cluster_nums.iter().copied());
        let mut pointers = vec![0; BRANCH_NUM];
        let mut merged_asn = Vec::with_capacity(seqs.len());
        let mut merged_post = Vec::with_capacity(seqs.len());
        for (asn, ps) in asn.into_iter().zip(pss) {
            let in_asn = recurred_clusterings[asn].0[pointers[asn]];
            let in_ps = &recurred_clusterings[asn].1[pointers[asn]];
            pointers[asn] += 1;
            merged_asn.push(offsets[asn] + in_asn);
            let mut posterior = Vec::with_capacity(total_cluster_num);
            for (p, &num) in ps.iter().zip(cluster_nums.iter()) {
                let lk = p - (num as f64).ln();
                posterior.extend(std::iter::repeat(lk).take(num));
            }
            for (t, p) in in_ps.iter().enumerate() {
                posterior[t + offsets[asn]] += p + (cluster_nums[asn] as f64).ln();
            }
            let sum: f64 = posterior.iter().map(|x| x.exp()).sum();
            assert!((1f64 - sum).abs() < 0.0001);
            merged_post.push(posterior);
        }
        (merged_asn, merged_post, total_score, total_cluster_num)
    }
}

fn cumsum<I: std::iter::Iterator<Item = usize>>(iter: I) -> Vec<usize> {
    iter.scan(0, |sum, x| {
        *sum += x;
        Some(*sum - x)
    })
    .collect()
}

fn filter_sub_clusters<'a>(
    seqs: &[&'a [u8]],
    ops: &[Vec<kiley::Op>],
    strands: &[bool],
    asn: &[usize],
    cluster: usize,
) -> (Vec<&'a [u8]>, Vec<Vec<kiley::Op>>, Vec<bool>) {
    let num: usize = asn.iter().filter(|&&x| x == cluster).count();
    let (mut sub_seqs, mut sub_ops, mut sub_strands) = (
        Vec::with_capacity(num),
        Vec::with_capacity(num),
        Vec::with_capacity(num),
    );
    for (((seq, ops), &strand), &asn) in seqs.iter().zip(ops).zip(strands).zip(asn.iter()) {
        if asn == cluster {
            sub_seqs.push(*seq);
            sub_ops.push(ops.clone());
            sub_strands.push(strand);
        }
    }
    (sub_seqs, sub_ops, sub_strands)
}

fn estim_copy_num(asn: &[usize], k: usize, copy_num: usize, coverage: f64) -> Vec<usize> {
    assert!(k <= copy_num, "{},{}", k, copy_num);
    let mut counts = vec![0f64; k];
    for &x in asn.iter() {
        counts[x] += 1f64;
    }
    let mut copy_numbers = vec![1; k];
    for _ in k..copy_num {
        *copy_numbers
            .iter_mut()
            .zip(counts.iter())
            .map(|(cp, cov)| ((cov - coverage * *cp as f64).powi(2), cp))
            .max_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
            .unwrap()
            .1 += 1;
    }
    let sum: usize = copy_numbers.iter().sum();
    assert_eq!(sum, copy_num);
    copy_numbers
}

fn update_by_clusterings(
    units: &mut [&mut Node],
    asn: &[usize],
    ops: &[Vec<kiley::Op>],
    pss: &[Vec<f64>],
) {
    for (node, ps) in units.iter_mut().zip(pss) {
        node.posterior.clear();
        node.posterior.extend(ps);
    }
    for (node, &asn) in units.iter_mut().zip(asn) {
        node.cluster = asn as u64;
    }
    for (node, ops) in units.iter_mut().zip(ops) {
        node.cigar = crate::misc::kiley_op_to_ops(ops);
    }
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

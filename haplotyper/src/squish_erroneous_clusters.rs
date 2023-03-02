use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct SquishConfig {
    ari_thr: f64,
    match_score: f64,
    mismatch_score: f64,
    count_thr: usize,
}
const BIAS_THR: f64 = 0.2;
impl SquishConfig {
    pub const fn new(
        ari_thr: f64,
        count_thr: usize,
        match_score: f64,
        mismatch_score: f64,
    ) -> Self {
        Self {
            ari_thr,
            count_thr,
            match_score,
            mismatch_score,
        }
    }
}

impl std::default::Default for SquishConfig {
    fn default() -> Self {
        Self {
            ari_thr: 0.5,
            match_score: 4f64,
            mismatch_score: -1f64,
            count_thr: 10,
        }
    }
}

pub trait SquishErroneousClusters {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig);
}

impl SquishErroneousClusters for DataSet {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig) {
        let chunk_class = classify_chunks(self, config);
        self.selected_chunks
            .iter_mut()
            .filter(|n| matches!(chunk_class.get(&n.id), Some(RelClass::Suspicious)))
            .for_each(|n| n.cluster_num = 1);
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .filter(|n| matches!(chunk_class.get(&n.chunk), Some(RelClass::Suspicious)))
            .for_each(|n| {
                n.cluster = 0;
                n.posterior = vec![0f64];
            });
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RelClass {
    Stiff,
    Isolated,
    Suspicious,
}

impl std::fmt::Display for RelClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RelClass::Stiff => write!(f, "Stiff"),
            RelClass::Isolated => write!(f, "Isolated"),
            RelClass::Suspicious => write!(f, "Suspicious"),
        }
    }
}

fn classify_chunks(ds: &DataSet, config: &SquishConfig) -> HashMap<u64, RelClass> {
    let mut chunk_pairs: HashMap<_, usize> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let nodes = read.nodes.iter().enumerate();
        for (i, n1) in nodes.filter(|n| n.1.is_biased(BIAS_THR)) {
            let n2s = read.nodes.iter().skip(i + 1);
            for n2 in n2s.filter(|n| n.is_biased(BIAS_THR)) {
                let key = (n1.chunk.min(n2.chunk), n1.chunk.max(n2.chunk));
                *chunk_pairs.entry(key).or_default() += 1;
            }
        }
    }
    let chunks: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|n| (n.id, n.cluster_num))
        .collect();
    chunk_pairs.retain(|_, val| config.count_thr < *val);
    chunk_pairs.retain(|(u1, u2), _| 1 < chunks[u1] && 1 < chunks[u2]);
    let adj_rand_indices: Vec<_> = chunk_pairs
        .par_iter()
        .map(|(&(u1, u2), _)| {
            let (cl1, cl2) = (chunks[&u1], chunks[&u2]);
            let (rel, count) = check_correl(ds, (u1, cl1), (u2, cl2));
            (u1, u2, (rel, count))
        })
        .collect();
    let mut touch_chunks: HashMap<_, Vec<_>> = HashMap::new();
    for (&(u1, u2), _) in chunk_pairs.iter() {
        touch_chunks.entry(u1).or_default().push(u2);
    }
    let stiff_chunks = match adj_rand_indices.is_empty() {
        true => HashSet::new(),
        false => classify(&adj_rand_indices, config),
    };
    // let stiff_chunks = classify_chunks_by_graph(ds, config);
    // ds.selected_chunks
    //     .iter()
    //     .map(|c| {
    //         let (is_stiff, touch_stiff) = match stiff_chunks.get(&c.id) {
    //             Some((is_stiff, touch_chunks)) => {
    //                 let touch_stiff = touch_chunks.iter().any(|c| match stiff_chunks.get(c) {
    //                     Some(x) => x.0,
    //                     None => false,
    //                 });
    //                 (*is_stiff, touch_stiff)
    //             }
    //             None => (false, false),
    //         };
    //         if is_stiff || 2 < c.copy_num {
    //             (c.id, RelClass::Stiff)
    //         } else if touch_stiff {
    //             (c.id, RelClass::Suspicious)
    //         } else {
    //             (c.id, RelClass::Isolated)
    //         }
    //     })
    //     .collect()
    ds.selected_chunks
        .iter()
        .map(|c| {
            let touch_stiff = touch_chunks.get(&c.id);
            let touch_stiff =
                touch_stiff.map(|rels| rels.iter().any(|to| stiff_chunks.contains(to)));
            if stiff_chunks.contains(&c.id) || 2 < c.copy_num {
                (c.id, RelClass::Stiff)
            } else if touch_stiff == Some(true) {
                let (to, ari, count) = adj_rand_indices
                    .iter()
                    .filter_map(|&(from, to, (ari, count))| {
                        if c.id == from {
                            Some((to, ari, count))
                        } else if c.id == to {
                            Some((from, ari, count))
                        } else {
                            None
                        }
                    })
                    .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
                    .unwrap();
                debug!("SUSPIC\t{}\t{to}\t{ari:.3}\t{count}", c.id);
                (c.id, RelClass::Suspicious)
            } else {
                (c.id, RelClass::Isolated)
            }
        })
        .collect()
}

// fn classify_chunks_by_graph(ds: &DataSet, config: &SquishConfig) -> HashMap<u64, (bool, Vec<u64>)> {
//     let mut chunk_pairs: HashMap<_, usize> = HashMap::new();
//     for read in ds.encoded_reads.iter() {
//         let nodes = read.nodes.iter().enumerate();
//         for (i, n1) in nodes.filter(|n| n.1.is_biased(BIAS_THR)) {
//             let n2s = read.nodes.iter().skip(i + 1);
//             for n2 in n2s.filter(|n| n.is_biased(BIAS_THR)) {
//                 let key = (n1.chunk.min(n2.chunk), n1.chunk.max(n2.chunk));
//                 *chunk_pairs.entry(key).or_default() += 1;
//             }
//         }
//     }
//     let chunks: HashMap<_, _> = ds
//         .selected_chunks
//         .iter()
//         .map(|n| (n.id, n.cluster_num))
//         .collect();
//     chunk_pairs.retain(|_, val| config.count_thr < *val);
//     chunk_pairs.retain(|(u1, u2), _| 1 < chunks[u1] && 1 < chunks[u2]);
//     let adj_rand_indices: Vec<_> = chunk_pairs
//         .par_iter()
//         .map(|(&(u1, u2), _)| {
//             let (cl1, cl2) = (chunks[&u1], chunks[&u2]);
//             let (rel, count) = check_correl(ds, (u1, cl1), (u2, cl2));
//             (u1, u2, (rel, count))
//         })
//         .collect();
//     let mut touch_chunks: HashMap<_, _> =
//         ds.selected_chunks.iter().map(|c| (c.id, vec![])).collect();
//     for (&(u1, u2), _) in chunk_pairs.iter() {
//         touch_chunks.entry(u1).or_default().push(u2);
//     }
//     let stiff_chunks = match adj_rand_indices.is_empty() {
//         true => HashSet::new(),
//         false => classify(&adj_rand_indices, config),
//     };
//     touch_chunks
//         .into_iter()
//         .map(|(c, chunks)| {
//             let is_stiff = stiff_chunks.contains(&c);
//             (c, (is_stiff, chunks))
//         })
//         .collect()
// }

fn check_correl(
    ds: &DataSet,
    (chunk1, cl1): (u64, usize),
    (chunk2, cl2): (u64, usize),
) -> (f64, usize) {
    let (mut c1, mut c2) = (vec![], vec![]);
    for read in ds.encoded_reads.iter() {
        let node1 = read
            .nodes
            .iter()
            .filter(|n| n.chunk == chunk1 && n.is_biased(BIAS_THR))
            .map(|n| n.cluster as usize)
            .min();
        let node2 = read
            .nodes
            .iter()
            .filter(|n| n.chunk == chunk2 && n.is_biased(BIAS_THR))
            .map(|n| n.cluster as usize)
            .min();
        if let (Some(n1), Some(n2)) = (node1, node2) {
            c1.push(n1);
            c2.push(n2);
        }
    }
    if c1.is_empty() {
        return (0f64, c1.len());
    }
    let c1_is_same = c1.iter().all(|&x| x == c1[0]);
    let c2_is_same = c2.iter().all(|&x| x == c2[0]);
    let rel_value = match (c1_is_same && c2_is_same, cl1 == 1 && cl2 == 1) {
        (true, true) => 0f64,
        (true, false) => 1f64,
        (false, _) => crate::misc::adjusted_rand_index(&c1, &c2),
    };
    if rel_value.is_nan() {
        warn!("\n{:?}\n{:?}", c1, c2);
        return (0f64, c1.len());
    }
    (rel_value, c1.len())
}

fn classify(adj_rand_indices: &[(u64, u64, (f64, usize))], config: &SquishConfig) -> HashSet<u64> {
    let mut nodes: HashMap<_, usize> = HashMap::new();
    for &(from, to, _) in adj_rand_indices.iter() {
        let len = nodes.len();
        nodes.entry(from).or_insert(len);
        let len = nodes.len();
        nodes.entry(to).or_insert(len);
    }
    let mut graph = vec![vec![]; nodes.len()];
    for &(from, to, (ari, count)) in adj_rand_indices.iter() {
        let ari = ari.max(0f64).min(1f64);
        let (from, to) = (nodes[&from], nodes[&to]);
        graph[from].push((to, ari, count));
        graph[to].push((from, ari, count));
    }
    let param = ClassifyParam::new(config.ari_thr, config.mismatch_score, config.match_score);
    let assignments = classify_nodes(&graph, nodes.len(), &param);
    nodes
        .iter()
        .filter_map(|(&uid, &node)| assignments[node].then_some(uid))
        .collect()
}

type RelGraph = Vec<Vec<(usize, f64, usize)>>;

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
fn classify_nodes(graph: &RelGraph, nodes: usize, param: &ClassifyParam) -> Vec<bool> {
    let mut assignments = vec![true; nodes];
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(3093240);
    for _t in 0..10 {
        wipe_through(graph, &mut assignments, param);
        mcmc(graph, &mut assignments, param, &mut rng);
    }
    assignments
}

#[derive(Debug, Clone)]
struct ClassifyParam {
    switch_point: f64,
    error_score: f64,
    correct_score: f64,
}

impl ClassifyParam {
    fn score(&self, ari: f64, count: usize) -> f64 {
        match ari <= self.switch_point {
            true => self.error_score * count as f64,
            false => self.correct_score * count as f64,
        }
    }
    fn new(switch_point: f64, error_score: f64, correct_score: f64) -> Self {
        Self {
            switch_point,
            error_score,
            correct_score,
        }
    }
}

impl std::fmt::Display for ClassifyParam {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:.3}\t{:.3}\t{:.3}",
            self.switch_point, self.error_score, self.correct_score
        )
    }
}

fn wipe_through(graph: &RelGraph, assignments: &mut [bool], param: &ClassifyParam) {
    let len = assignments.len();
    let prev_score = assignment_score(graph, assignments, param);
    for i in 0..len {
        if 0f64 < diff_on_flip(graph, assignments, i, param) {
            assignments[i] = !assignments[i];
        }
    }
    let after_score = assignment_score(graph, assignments, param);
    trace!("WIPE\t{prev_score:.3}\t{after_score:.3}")
}

fn diff_on_flip(
    graph: &RelGraph,
    assignments: &[bool],
    target: usize,
    param: &ClassifyParam,
) -> f64 {
    let sum_weight = graph[target]
        .iter()
        .filter(|&&(to, _, _)| assignments[to])
        .map(|&(_, ari, count)| param.score(ari, count))
        .sum::<f64>();
    match assignments[target] {
        true => -sum_weight,
        false => sum_weight,
    }
}
use rand::Rng;
fn mcmc<R: Rng>(graph: &RelGraph, assignments: &mut [bool], param: &ClassifyParam, rng: &mut R) {
    let prev_score = assignment_score(graph, assignments, param);
    for _ in 0..1000 {
        let i = rng.gen_range(0..assignments.len());
        let diff = diff_on_flip(graph, assignments, i, param);
        let prob = diff.min(0f64).exp();
        if rng.gen_bool(prob) {
            assignments[i] = !assignments[i];
        }
    }
    let after_score = assignment_score(graph, assignments, param);
    trace!("MCMC\t{prev_score:.3}\t{after_score:.3}");
}

fn assignment_score(graph: &RelGraph, assignments: &[bool], param: &ClassifyParam) -> f64 {
    graph
        .iter()
        .enumerate()
        .map(|(from, edges)| {
            edges
                .iter()
                .filter(|&&(to, _, _)| assignments[to] && assignments[from] && from < to)
                .map(|&(_, ari, count)| param.score(ari, count))
                .sum::<f64>()
        })
        .sum()
}

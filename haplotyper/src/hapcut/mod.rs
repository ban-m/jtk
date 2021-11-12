//! A module to implement the original HapCUT algorithm (MEC-based) in a chunk-grpah.
use log::*;
use rand::{prelude::SliceRandom, Rng, SeedableRng};
use rand_distr::num_traits::Signed;
use rand_xoshiro::Xoshiro256Plus;

use definitions::*;
use std::collections::{HashMap, HashSet};
pub trait HapCut {
    /// Run HapCUT algorithm. Return the consistency score of the result.
    fn hapcut(&mut self, c: &HapCutConfig) -> (i64, Vec<Phase>);
    fn hapcut_squish_diplotig(&mut self, c: &HapCutConfig) -> i64;
    /// Run HapCUT2 algorithm. Return the (scaled-)likelihood of the result.
    fn hapcut2(&mut self, c: &HapCut2Config) -> f64;
}

#[derive(Debug, Clone)]
pub struct HapCutConfig {
    seed: u64,
    // // Reward for match. Positive.
    // match_reward: i64,
    // // Mismatch penalty. Negative.
    // mism_penalty: i64,
    // // Penalty for dont care position.
    // dont_care_score: i64,
    // Rounding threshold.
    // For k-dim weight vector x, if x[i].abs() < round_thr holds,
    // the corresponding phase would be ``graph'', or, regardedas diplochunk.
    // Usually, this threshold is less than 1/4.
    round_threshold: f64,
    // The count threshold used to construct the adj graph.
    // The edges with count less than or equal to this value would be dicarded.
    count_thr: u32,
    // Copy number information
    copy_number: HashMap<u64, usize>,
}

impl HapCutConfig {
    pub fn set_thr(&mut self, thr: f64) {
        self.round_threshold = thr;
    }
    pub fn set_seed(&mut self, seed: u64) {
        self.seed = seed;
    }
    fn update_seed(&mut self) {
        self.seed += 1;
    }
    fn seed(&self) -> u64 {
        self.seed
    }
    pub fn new(
        seed: u64,
        round_threshold: f64,
        count_thr: u32,
        copy_number: HashMap<u64, usize>,
    ) -> Self {
        Self {
            seed,
            round_threshold,
            count_thr,
            copy_number,
        }
    }
}

impl std::default::Default for HapCutConfig {
    fn default() -> Self {
        Self {
            seed: 0,
            // match_reward: 1,
            // mism_penalty: -1,
            // dont_care_score: 0,
            round_threshold: 0.25,
            count_thr: 0,
            copy_number: HashMap::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct HapCut2Config {}

#[derive(Clone)]
pub struct Phase {
    unit: u64,
    cluster_size: usize,
    // -1,0,1
    phase: Vec<i8>,
}

impl std::fmt::Debug for Phase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{:?}", self.unit, self.phase)
    }
}

impl std::ops::Index<u64> for Phase {
    type Output = i8;
    fn index(&self, index: u64) -> &Self::Output {
        self.phase.index(index as usize)
    }
}

impl Phase {
    pub fn phases(&self) -> Vec<((u64, u64), i8)> {
        self.phase
            .iter()
            .enumerate()
            .map(|(cl, phase)| ((self.unit, cl as u64), *phase))
            .collect()
    }
    pub fn phase(&self) -> &[i8] {
        &self.phase
    }
    pub fn new(unit: u64, phase: Vec<i8>) -> Self {
        Self {
            cluster_size: phase.len(),
            unit,
            phase,
        }
    }
    fn flip_cut(phases: &mut [Phase], cut: &[bool]) {
        phases
            .iter_mut()
            .zip(cut.iter())
            .filter(|&(_, &c)| c)
            .for_each(|(p, _)| p.flip());
    }
    fn flip(&mut self) {
        self.phase.iter_mut().for_each(|x| *x *= -1);
    }
    fn new_with_rand<R: Rng>(rng: &mut R, unit: u64, cluster_size: usize) -> Self {
        //         const INIT_WEIGHT: [f64; 3] = [0.45, 0.1, 0.45];
        // const CHOISES: [usize; 3] = [0, 1, 2];
        // let phase: Vec<_> = (0..cluster_size)
        //     .map(|_| [-1, 0, 1][*CHOISES.choose_weighted(rng, |&k| INIT_WEIGHT[k]).unwrap()])
        //     .collect();
        let mut phase: Vec<_> = (0..cluster_size)
            .map(|i| if i < cluster_size / 2 { -1 } else { 1 })
            .collect();
        phase.shuffle(rng);
        // if cluster_size == 1 {
        if cluster_size != 2 {
            phase = vec![0; cluster_size];
        }
        Self {
            unit,
            cluster_size,
            phase,
        }
    }
    // Return how much consistencies would be changed on fliping haplotyping on chunks.
    // i->the consistency change by flipping the i-th phase.
    pub fn consistency_change_by_flip(
        phases: &[Phase],
        reads: &[(u64, Vec<(u64, u64)>)],
        c: &HapCutConfig,
    ) -> Vec<i64> {
        let mut consistency_changes = vec![0; phases.len()];
        for (_, read) in reads.iter() {
            let (plus, minus) = Self::get_consis_on_phases(phases, read, c);
            let consis = plus.min(minus);
            for &(unit, cl) in read {
                // If we flip the phase on phases[unit]...
                consistency_changes[unit as usize] += match phases[unit as usize][cl] {
                    1 => (plus + 2).min(minus - 2) - consis,
                    -1 => (plus - 2).min(minus + 2) - consis,
                    _ => 0,
                };
            }
        }
        consistency_changes.iter_mut().for_each(|x| *x /= 2);
        consistency_changes
    }
    // Return consistency score of the current phase.
    // *Smaller is better*.
    fn consistency(phases: &[Phase], reads: &[(u64, Vec<(u64, u64)>)], c: &HapCutConfig) -> i64 {
        let total_consistency: i64 = reads
            .iter()
            .map(|(_, read)| Phase::consistency_of_read(phases, read, c))
            .sum();
        total_consistency / 2
    }
    fn consistency_of_read(phases: &[Phase], read: &[(u64, u64)], c: &HapCutConfig) -> i64 {
        let (plus, minus) = Self::get_consis_on_phases(phases, read, c);
        plus.min(minus)
    }
    fn get_consis_on_phases(
        phases: &[Phase],
        read: &[(u64, u64)],
        _c: &HapCutConfig,
    ) -> (i64, i64) {
        // Match -> 0, Mismatch -> 2, Diplochunk -> 1
        let plus: i64 = read
            .iter()
            .map(|&(unit, cl)| 1 - phases[unit as usize][cl] as i64)
            .sum();
        // Match -> 2, Mismatch -> 0, Diplochunk -> 1 (Corresponding to the other side of the haplotype).
        let minus: i64 = read
            .iter()
            .map(|&(unit, cl)| 1 + phases[unit as usize][cl] as i64)
            .sum();
        (plus, minus)
    }
    // Optimizing a unit.
    // It re-assign the haplotype assignment on each chunk, by choosing a ``near-parallel'' integer vector
    // with `cluster_size`-dimensional vector `w`, where the k-th element of `w` is defined as
    // the sum of the w_rk over the all reads having (unit,k) in them.
    // w_rk :=
    // (sum of phases[u][c] of the read r)/(length of r - 1) divided by
    // (sum of abs(phases[u][c]) of the read r)/(length of r - 1)
    #[allow(dead_code)]
    fn optimize_unit(phases: &mut [Phase], reads: &[(u64, Vec<(u64, u64)>)], c: &HapCutConfig) {
        let weight_vector: Vec<_> = phases
            .iter()
            .map(|phase| {
                let (mut ws, aws) = Self::get_weight_vector(phase, phases, reads);
                ws.iter_mut().zip(aws).for_each(|(w, a)| *w /= a);
                ws
            })
            .collect();
        phases
            .iter_mut()
            .zip(weight_vector.iter())
            .for_each(|(phase, weights)| {
                assert_eq!(weights.len(), phase.phase.len());
                phase.phase = Self::approx_parallel_vector(weights, c);
            });
    }
    // Same as `optimize_unit`, but it only re-assign phases on the multi-copy(>2) chunks.
    fn optimize_unit_mc(phases: &mut [Phase], reads: &[(u64, Vec<(u64, u64)>)], c: &HapCutConfig) {
        let weight_vector: Vec<_> = phases
            .iter()
            .filter(|p| 2 < p.cluster_size)
            .map(|phase| {
                let (mut ws, aws) = Self::get_weight_vector(phase, phases, reads);
                ws.iter_mut().zip(aws).for_each(|(w, a)| *w /= a);
                ws
            })
            .collect();
        phases
            .iter_mut()
            .filter(|p| 2 < p.cluster_size)
            .zip(weight_vector.iter())
            .for_each(|(phase, weights)| {
                assert_eq!(weights.len(), phase.phase.len());
                phase.phase = Self::approx_parallel_vector(weights, c);
            });
    }
    // TODO: This function can be faster, linear search on read in not neede.(?)
    fn get_weight_vector(
        phase: &Phase,
        phases: &[Phase],
        reads: &[(u64, Vec<(u64, u64)>)],
    ) -> (Vec<f64>, Vec<f64>) {
        let mut weight_vector: Vec<_> = vec![0f64; phase.cluster_size];
        let mut weights_abs: Vec<_> = vec![0f64; phase.cluster_size];
        for (_, read) in reads.iter() {
            for node in read.iter().filter(|&&(n, _)| n == phase.unit) {
                let read = read
                    .iter()
                    .filter(|&n| n != node)
                    .map(|&(unit, cluster)| phases[unit as usize][cluster] as i32);
                let total: i32 = read.clone().sum();
                let abs_total: i32 = read.clone().map(|x| x.abs()).sum();
                let num = read.count();
                if num != 0 {
                    weight_vector[node.1 as usize] += total as f64 / num as f64;
                    weights_abs[node.1 as usize] += abs_total as f64 / num as f64;
                }
            }
        }
        (weight_vector, weights_abs)
    }
    // Given vector xs, return {-1,0,1}^k vector ys maximizing <xs,ys>.
    // Normally, it is 1 if xs[k] is positive, -1 otherwise.
    // However, to polish the signature, we map small xs[k] value into 0.
    // This threshold can be tune by c.
    fn approx_parallel_vector(weights: &[f64], c: &HapCutConfig) -> Vec<i8> {
        weights
            .iter()
            .map(|x| match c.round_threshold < x.abs() {
                true if x.is_sign_positive() => 1,
                true => -1,
                false => 0,
            })
            .collect()
    }
    // Optimize the phasing by finding a positive score cut.
    // If the cut is found and improved, return true, otherwise, false.
    fn optimize_by_cut(
        phases: &mut [Phase],
        reads: &[(u64, Vec<(u64, u64)>)],
        c: &HapCutConfig,
    ) -> (Vec<bool>, f64) {
        let consistency_graph = Self::to_consistency_graph(phases, reads, c);
        // for (i, edges) in consistency_graph.iter().enumerate() {
        //     if edges.iter().any(|e| e.1 > 0f64) {
        //         eprintln!("{}\t{:?}", i, edges);
        //     }
        // }
        // eprintln!("============");
        Self::find_cut(&consistency_graph, &c)
    }
    // Make the adjacency graph for the given dataset.
    fn to_consistency_graph(
        phases: &[Phase],
        reads: &[(u64, Vec<(u64, u64)>)],
        c: &HapCutConfig,
    ) -> Vec<Vec<(usize, f64)>> {
        // i -> vector of (read id, Hi[ri], |r|)
        let pileups = Self::to_pileups(phases, reads, c);
        let connections = Self::get_connections(phases.len(), reads, c);
        connections
            .iter()
            .enumerate()
            .map(|(from, edges)| {
                edges
                    .iter()
                    .map(|&to| (to, Self::get_edge_weight(&pileups[from], &pileups[to])))
                    .collect()
            })
            .collect()
    }
    // This is |p1| * |p2| algorithm.
    // We can make it O(|p1| + |p2|) by using peekable iterator,
    // but it can be tricky under the situation when there is multiple occurence of the
    // same unit, with different cluster, i.e., in a loop region.
    fn get_edge_weight(p1: &[(u64, i8, usize)], p2: &[(u64, i8, usize)]) -> f64 {
        p1.iter()
            .filter(|&&(_, _, len)| 1 < len)
            .map(|&(id1, phase1, len)| {
                let weight: i8 = p2
                    .iter()
                    .filter(|&&(id2, _, _)| id1 == id2)
                    .map(|(_, phase2, _)| phase1 * phase2)
                    .sum();
                -weight as f64 / (len - 1) as f64
            })
            .sum()
    }
    // i -> vector of (read id, Hi[ri], |r|)
    fn to_pileups(
        phases: &[Phase],
        reads: &[(u64, Vec<(u64, u64)>)],
        _c: &HapCutConfig,
    ) -> Vec<Vec<(u64, i8, usize)>> {
        let mut pileups: Vec<_> = vec![vec![]; phases.len()];
        for &(id, ref read) in reads.iter() {
            let len = read.len();
            for (unit, cluster) in read.iter().map(|&(u, c)| (u as usize, c)) {
                pileups[unit].push((id, phases[unit][cluster], len));
            }
        }
        pileups
            .iter_mut()
            .for_each(|pileup| pileup.sort_by_key(|x| x.0));
        pileups
    }
    // Get a simple adjacency list.
    fn get_connections(
        len: usize,
        reads: &[(u64, Vec<(u64, u64)>)],
        c: &HapCutConfig,
    ) -> Vec<Vec<usize>> {
        let mut edges_and_counts: Vec<HashMap<usize, u32>> = vec![HashMap::new(); len];
        for (_, read) in reads.iter() {
            for (i, &(from, _)) in read.iter().enumerate() {
                let from = from as usize;
                for &(to, _) in read.iter().skip(i + 1) {
                    let to = to as usize;
                    *edges_and_counts[from].entry(to).or_default() += 1;
                    *edges_and_counts[to].entry(from).or_default() += 1;
                }
            }
        }
        edges_and_counts
            .iter()
            .map(|counts| {
                let mut edges: Vec<_> = counts
                    .iter()
                    .filter_map(|(&node, &count)| (c.count_thr < count).then(|| node))
                    .collect();
                edges.sort_unstable();
                edges
            })
            .collect()
    }
    // Find a cut, return cut and its capacity.
    fn find_cut(graph: &[Vec<(usize, f64)>], c: &HapCutConfig) -> (Vec<bool>, f64) {
        // Serialized version of the graph.
        let mut edges = Vec::with_capacity(graph.iter().fold(0, |x, es| x + es.len()));
        for (from, es) in graph.iter().enumerate() {
            edges.extend(es.iter().map(|&(to, w)| (from, to, w)));
        }
        let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(c.seed);
        let trial_num = 10;
        let init_edges: Vec<(usize, usize, f64)> = (0..trial_num)
            .filter_map(|_| {
                edges
                    .choose_weighted(&mut rng, |(_, _, w)| w.max(0f64))
                    .ok()
            })
            .copied()
            .collect();
        let seed = c.seed();
        init_edges
            .iter()
            .enumerate()
            .map(|(i, &edge)| Phase::propose_cut(graph, edge, c, seed + i as u64 + 1))
            .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
            .unwrap_or_else(|| (vec![], 0f64))
    }
    // Return A(v) = w(v,S_1)/|e(v,S_1)| - w(v,S_2)/|e(v,S_2)|.
    // Intuitively, if A(v) > 0, v should be in S_2, and S_1 otherwise.
    fn node_capacity(node: usize, graph: &[Vec<(usize, f64)>], cut: &[i8]) -> f64 {
        let mut weight = 0f64;
        for &(to, w) in graph[node].iter() {
            if cut[to] == 1 {
                weight += w;
            } else if cut[to] == -1 {
                weight -= w;
            }
        }
        weight
    }
    // Gen -1 or 1 evenly.
    fn gen<R: Rng>(rng: &mut R) -> i8 {
        match rng.gen_bool(0.5) {
            true => 1,
            false => -1,
        }
    }
    // return w(S_1,S_2).
    fn cut_capacity(graph: &[Vec<(usize, f64)>], cut: &[i8]) -> f64 {
        graph
            .iter()
            .zip(cut.iter())
            .map(|(edges, &belong)| {
                edges
                    .iter()
                    .filter(|&&(to, _)| belong != cut[to])
                    .fold(0f64, |x, &(_, w)| x + w)
            })
            .sum()
    }
    fn propose_cut(
        graph: &[Vec<(usize, f64)>],
        (from, to, _): (usize, usize, f64),
        _c: &HapCutConfig,
        seed: u64,
    ) -> (Vec<bool>, f64) {
        let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(seed);
        // 1->S_1, -1->S_2, 0->not searced yet.
        let mut cut: Vec<i8> = vec![0; graph.len()];
        // eprintln!("CUT:{}-{}", from, to);
        cut[from] = 1;
        cut[to] = -1;
        loop {
            let next_target = cut
                .iter()
                .enumerate()
                .filter(|&(_, &belonging)| belonging == 0)
                .map(|(node, _)| (node, Self::node_capacity(node, graph, &cut)))
                .max_by(|x, y| (x.1.abs()).partial_cmp(&y.1.abs()).unwrap());
            // if let Some((n, w)) = next_target.clone() {
            //     eprintln!("PUSH\t{}\t{}", n, w);
            // }
            match next_target {
                Some((node, weight)) if weight < 0f64 => cut[node] = 1,
                Some((node, weight)) if weight > 0f64 => cut[node] = -1,
                Some((node, _)) => cut[node] = Self::gen(&mut rng),
                None => break,
            }
        }
        // let (cut1, cut2): (Vec<_>, Vec<_>) = (0..cut.len()).partition(|&x| cut[x] == 1);
        // eprintln!("CUT1\t{:?}", cut1);
        // eprintln!("CUT2\t{:?}", cut2);
        // Fine tune the cut.
        let mut cut_capacity = Self::cut_capacity(graph, &cut);
        // TODO: Is this OK? Should be continue this loop
        // each time after we change the cut,
        // rather than chaning the cut in a batch.
        loop {
            // Preverse the previous cut.
            let prev_cut = cut.clone();
            for (node, belong) in cut.iter_mut().enumerate() {
                *belong = match Self::node_capacity(node, graph, &prev_cut) {
                    w if w < 0f64 => 1,
                    w if w > 0f64 => -1,
                    _ => Self::gen(&mut rng),
                };
            }
            let new_capacity = Self::cut_capacity(graph, &cut);
            if cut_capacity < new_capacity {
                cut_capacity = new_capacity;
            } else {
                // Revert & break.
                cut = prev_cut;
                break;
            }
        }
        assert!(cut.iter().all(|&x| x == 1 || x == -1));
        let cut: Vec<_> = cut.iter().map(|&x| x == 1).collect();
        (cut, cut_capacity)
    }
}

impl HapCut for DataSet {
    fn hapcut_squish_diplotig(&mut self, c: &HapCutConfig) -> i64 {
        let reads: Vec<(u64, Vec<(u64, u64)>)> = self
            .encoded_reads
            .iter()
            .map(|r| (r.id, r.nodes.iter().map(|x| (x.unit, x.cluster)).collect()))
            .collect();
        let (phases, consistency) = hapcut(&reads, c);
        // Squishing nodes.
        let is_uninformative: HashSet<_> = phases
            .iter()
            .enumerate()
            .filter_map(|(i, phase)| {
                assert_eq!(phase.unit, i as u64);
                (phase.phase.iter().all(|&p| p == 0)).then(|| i as u64)
            })
            .collect();
        debug!("SQUISHED\t{}\tCLUSTER", is_uninformative.len());
        self.encoded_reads
            .iter_mut()
            .flat_map(|read| read.nodes.iter_mut())
            .filter(|node| is_uninformative.contains(&node.unit))
            .for_each(|node| node.cluster = 0);
        // self.encoded_reads
        //     .iter_mut()
        //     .flat_map(|read| read.nodes.iter_mut())
        //     .for_each(|node| {
        //         node.cluster = mapping[node.unit as usize][node.cluster as usize];
        //     });
        // This is a vector of the hashset.
        consistency
    }
    fn hapcut(&mut self, c: &HapCutConfig) -> (i64, Vec<Phase>) {
        let reads: Vec<(u64, Vec<(u64, u64)>)> = self
            .encoded_reads
            .iter()
            .map(|r| (r.id, r.nodes.iter().map(|x| (x.unit, x.cluster)).collect()))
            .collect();
        let (phases, consistency) = hapcut(&reads, c);
        let changes = Phase::consistency_change_by_flip(&phases, &reads, c);
        trace!("HAPCUT\tChange\tunit\tdiff");
        for (phase, change) in phases.iter().zip(changes.iter()) {
            trace!("HAPCUT\tChange\t{}\t{}", phase.unit, change);
        }
        for phase in phases.iter() {
            trace!("HAPCUT\tDUMP\t{:?}", phase);
        }
        let phased_blocks = split_phase_blocks(&phases, &reads, c);
        for pb in phased_blocks.iter() {
            let pb: Vec<_> = pb.iter().map(|(u, cl)| format!("{}-{}", u, cl)).collect();
            trace!("HAPCUT\tPB\t{}", pb.join("\t"));
        }
        self.assignments = reads_into_phased_block(&reads, &phased_blocks, c);
        (consistency, phases)
    }
    fn hapcut2(&mut self, _c: &HapCut2Config) -> f64 {
        todo!()
    }
}

fn hapcut(reads: &[(u64, Vec<(u64, u64)>)], c: &HapCutConfig) -> (Vec<Phase>, i64) {
    debug!("HAPCUT\tStart\t1");
    // First, associate cluster/phase randomly.
    // Vector of phases. [i] -> i-th chunk.
    let mut phases = gen_random_phase(reads, c).unwrap();
    let mut c = c.clone();
    for t in 0..50 {
        trace!("HAPCUT\tIter\t{}\t1", t);
        let (cut, _) = Phase::optimize_by_cut(&mut phases, reads, &c);
        let prev_consis = Phase::consistency(&phases, reads, &c);
        trace!("HAPCUT\tCONSIS\t{}\t{}\t0", t, prev_consis);
        Phase::flip_cut(&mut phases, &cut);
        let consis = Phase::consistency(&phases, reads, &c);
        trace!("HAPCUT\tCONSIS\t{}\t{}\t1", t, consis);
        if prev_consis <= consis {
            Phase::flip_cut(&mut phases, &cut);
        }
        // if consis < prev_consis {
        // let cut: Vec<_> = cut
        //     .iter()
        //     .enumerate()
        //     .filter_map(|(i, &c)| c.then(|| i))
        //     .collect();
        // eprintln!("CUT\t{}\t{:?}", prev_consis - consis, cut);
        // }
        if 20 < t {
            Phase::optimize_unit_mc(&mut phases, reads, &c);
        }
        // if 220 < t {
        //     Phase::optimize_unit(&mut phases, reads, &c);
        // }
        // let consis = Phase::consistency(&phases, reads, &c);
        // trace!("HapCut\tCONSIS\t{}\t{}\t2", t, consis);
        c.update_seed();
    }
    // Phase::optimize_unit(&mut phases, reads, &c);
    let consis = Phase::consistency(&phases, reads, &c);
    (phases, consis)
}

fn get_number_ofchunk(reads: &[(u64, Vec<(u64, u64)>)]) -> Option<usize> {
    reads
        .iter()
        .filter_map(|(_, r)| r.iter().map(|x| x.0).max())
        .max()
        .map(|x| x as usize + 1)
}

fn gen_random_phase(reads: &[(u64, Vec<(u64, u64)>)], c: &HapCutConfig) -> Option<Vec<Phase>> {
    let chunk_num = get_number_ofchunk(reads)?;
    let mut max_cluster_num = vec![0; chunk_num];
    for (_, read) in reads.iter() {
        for (chunk, cluster) in read.iter().map(|&(n, c)| (n as usize, c as usize)) {
            if let Some(cl) = max_cluster_num.get_mut(chunk) {
                *cl = (*cl).max(cluster + 1);
            }
        }
    }
    let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(c.seed);
    let init_phases: Vec<_> = max_cluster_num
        .iter()
        .enumerate()
        .map(|(idx, &cl_num)| {
            let mut phase = Phase::new_with_rand(&mut rng, idx as u64, cl_num);
            if matches!(c.copy_number.get(&(idx as u64)), Some(1)) {
                // If this unit is haplotype specific, add phase information.
                assert_eq!(cl_num, 1);
                phase.phase = if rng.gen_bool(0.5) { vec![1] } else { vec![-1] };
            }
            phase
        })
        .collect();
    Some(init_phases)
}

/// Split phases into phased blocks. In other words, for each plus(1)/minus(-1) clusters,
/// compute connected components, and aggregate them as a vector of hashsets.
/// Note that all the gray(0) clusters would be removed after this splits.
/// They are regarded as ``untagged'' reads and should be rescued afterwords.
fn split_phase_blocks(
    phases: &[Phase],
    reads: &[(u64, Vec<(u64, u64)>)],
    c: &HapCutConfig,
) -> Vec<HashSet<(u64, u64)>> {
    let mut phase_blocks = vec![];
    for target in [-1, 1] {
        let mut hap_nodes = vec![];
        for phase in phases.iter() {
            for (cl, &hap) in phase.phase.iter().enumerate() {
                if hap == target {
                    hap_nodes.push((phase.unit, cl as u64));
                }
            }
        }
        phase_blocks.extend(connected_components(&hap_nodes, reads, c));
    }
    phase_blocks
}

fn connected_components(
    nodes: &[(u64, u64)],
    reads: &[(u64, Vec<(u64, u64)>)],
    _config: &HapCutConfig,
) -> Vec<HashSet<(u64, u64)>> {
    use crate::find_union::FindUnion;
    let nodes: HashMap<_, _> = nodes
        .iter()
        .copied()
        .enumerate()
        .map(|(x, y)| (y, x))
        .collect();
    let mut fu = FindUnion::new(nodes.len());
    for (_, read) in reads {
        let mut read = read.iter().filter(|n| nodes.contains_key(n));
        let mut from = match read.next().map(|n| nodes.get(n)) {
            Some(Some(idx)) => *idx,
            _ => continue,
        };
        for node in read {
            // This unwrap never panics.
            let to = *nodes.get(node).unwrap();
            fu.unite(from, to);
            from = to;
        }
    }
    let mut nodes: Vec<((u64, u64), usize)> = nodes.into_iter().collect();
    nodes.sort_by_key(|x| x.1);
    // O(KN)
    let mut phased_blocks = vec![];
    for i in 0..nodes.len() {
        if fu.find(i) == Some(i) {
            // i is the representative for a phased block.
            let phased_block: HashSet<_> = nodes
                .iter()
                .filter_map(|&(n, j)| (fu.find(j) == Some(i)).then(|| n))
                .collect();
            phased_blocks.push(phased_block);
        }
    }
    phased_blocks
}

/// Split reads into phased blocks. If a read does not contain any phased nodes,
/// Return None, otherwise, return the most probable phase block/assignment.
fn reads_into_phased_block(
    reads: &[(u64, Vec<(u64, u64)>)],
    phased_blocks: &[HashSet<(u64, u64)>],
    config: &HapCutConfig,
) -> Vec<Assignment> {
    reads
        .iter()
        .filter_map(|&(id, ref read)| {
            read_into_phased_block(read, &phased_blocks, config)
                .map(|block_id| Assignment::new(id, block_id))
        })
        .collect()
}

fn read_into_phased_block(
    read: &[(u64, u64)],
    phased_blocks: &[HashSet<(u64, u64)>],
    _config: &HapCutConfig,
) -> Option<usize> {
    let overlap = |hs: &HashSet<(u64, u64)>, read: &[(u64, u64)]| {
        read.iter().filter(|&elm| hs.contains(elm)).count()
    };
    phased_blocks
        .iter()
        .enumerate()
        .map(|(block_id, pb)| (block_id, overlap(pb, read)))
        .filter(|&(_, overlap)| 0 < overlap)
        .max_by_key(|x| x.1)
        .map(|x| x.0)
}

#[cfg(test)]
pub mod hapcut {
    use super::*;
    // Generate HapCut test dataset.
    // The first n fragments come from the hapA, others from hapB.
    // varnum: variant number,
    // len: length of the fragments, error: error rate, diprate: diploid probability.
    // Wrong clustering: fraction of uninformative variant.
    // loop_num: finally, we iteratively merge two variants as the `same` by`loop_num` times.
    fn gen_test_dataset(
        (num_reads, len, error): (usize, usize, f64),
        (varnum, diprate, wrong_cl, loop_num): (usize, f64, f64, usize),
        seed: u64,
    ) -> Vec<(u64, Vec<(u64, u64)>)> {
        let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(seed);
        // 0->usual, 1-> dip(cluster num = 1), 2-> error(cluster is uninformative).
        let is_dip_or_wrong: Vec<_> = (0..varnum)
            .map(|_| {
                if rng.gen_bool(diprate) {
                    1
                } else if rng.gen_bool(wrong_cl) {
                    2
                } else {
                    0
                }
            })
            .collect();
        eprintln!("{}", is_dip_or_wrong.len());
        let dip: Vec<_> = is_dip_or_wrong
            .iter()
            .enumerate()
            .filter_map(|(i, &ty)| (ty == 1).then(|| i))
            .collect();
        eprintln!("DIP\t{:?}", dip);
        let uninf: Vec<_> = is_dip_or_wrong
            .iter()
            .enumerate()
            .filter_map(|(i, &ty)| (ty == 2).then(|| i))
            .collect();
        eprintln!("Uninf\t{:?}", uninf);
        let mut gen = |i| -> Vec<(u64, u64)> {
            let hap = i < num_reads;
            let start = rng.gen_range(0..varnum - len + 1);
            let mut fragment: Vec<_> = is_dip_or_wrong
                .iter()
                .enumerate()
                .skip(start)
                .take(len)
                .map(|(var, is_d_w)| match is_d_w {
                    0 if rng.gen_bool(error) => (var as u64, (!hap) as u64),
                    0 => (var as u64, hap as u64),
                    1 => (var as u64, 0),
                    _ => (var as u64, rng.gen_bool(0.5) as u64),
                })
                .collect();
            if rng.gen_bool(0.5) {
                fragment.reverse();
            }
            fragment
        };
        let mut fragments: Vec<_> = (0..2 * num_reads).map(|i| (i as u64, gen(i))).collect();
        let mut merged: Vec<Vec<_>> = (0..varnum).map(|i| vec![i]).collect();
        for _ in 0..loop_num {
            let removed = merged.remove(rng.gen_range(0..merged.len()));
            let push_idx = rng.gen_range(0..merged.len());
            merged[push_idx].extend(removed);
        }
        for cl in merged.iter().filter(|xs| 1 < xs.len()) {
            eprintln!("Merge\t{:?}", cl);
        }
        // var idx -> offset of the variants.
        let mut map_to: HashMap<_, _> = HashMap::new();
        for vars in merged {
            let min = (*vars.iter().min().unwrap()) as u64;
            let mut haps_so_far = 0;
            for var in vars {
                map_to.insert(var as u64, (min, haps_so_far));
                haps_so_far += match is_dip_or_wrong[var] {
                    1 => 1, // Homologous chunk.
                    _ => 2, // Hap chunk.
                };
            }
        }
        for (_, frag) in fragments.iter_mut() {
            frag.iter_mut().for_each(|node| {
                let (new_var, offset) = map_to[&node.0];
                *node = (new_var, node.1 + offset);
            });
        }
        fragments
    }
    // read_params: (number, length, error rate)
    // var params: (number, diprate, uninformative rate, loop_number),
    // HapCut with length of 2. Error free.
    #[test]
    fn test1() {
        let config = HapCutConfig::default();
        let read_params = (100, 2, 0f64);
        let var_params = (10, 0f64, 0f64, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, consis) = hapcut(&reads, &config);
        eprintln!("{}", consis);
        assert_eq!(0, consis);
        let first_phase = &phases[2];
        for phase in phases.iter() {
            assert!(phase.phase == first_phase.phase);
        }
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let mut assignments: HashMap<_, Vec<_>> = HashMap::new();
        for Assignment { id, cluster } in reads_into_phased_block(&reads, &phased_blocks, &config) {
            assignments.entry(cluster).or_default().push(id);
        }
        assert_eq!(assignments.len(), 2);
        for members in assignments.values() {
            let answer = members[0] < read_params.0 as u64;
            for &mem in members.iter() {
                assert_eq!((mem < read_params.0 as u64), answer)
            }
        }
    }
    // HapCut with length of 2. Error = 0.1.
    #[test]
    fn test2() {
        let config = HapCutConfig::default();
        let read_params = (100, 2, 0.1);
        let var_params = (10, 0f64, 0f64, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, consis) = hapcut(&reads, &config);
        assert!(0 < consis);
        let first_phase = &phases[2];
        for phase in phases.iter() {
            assert!(phase.phase == first_phase.phase);
        }
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // HapCut with length of 3. Error = 0.1.
    #[test]
    fn test3() {
        let config = HapCutConfig::default();
        let read_params = (100, 3, 0.1);
        let var_params = (10, 0f64, 0f64, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let first_phase = &phases[2];
        for phase in phases.iter() {
            assert!(phase.phase == first_phase.phase);
        }
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // HapCut with length of 3. Error = 0.1, DipRate = 0.1
    #[test]
    fn test4() {
        let config = HapCutConfig::default();
        let read_params = (100, 3, 0.1);
        let var_params = (10, 0.1, 0f64, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            eprintln!("{:?}", phase);
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 8 / 10 < largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // HapCut with length of 3.
    // Error = 0.1, DipRate = 0.1, Uninformative = 0.1.
    #[test]
    fn test5() {
        let config = HapCutConfig::default();
        let read_params = (100, 3, 0.1);
        let var_params = (10, 0.1, 0.1, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 8 / 10 <= largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // Graph version of the each test
    // HapCut with length of 3.
    // Error = 0.1, DipRate = 0.1, Uninformative = 0.1, Two loops.
    #[test]
    fn test_graph() {
        let config = HapCutConfig::default();
        let read_params = (1000, 4, 0.1);
        let var_params = (50, 0.1, 0.1, 2);
        let seed = 38940;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            eprintln!("{:?}", phase);
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 7 / 10 <= largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        if phased_blocks.len() != 2 {
            for bl in phased_blocks.iter() {
                eprintln!("PB\t{:?}", bl);
            }
        }
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // Test large dataset.
    // HapCut with length of 10, error 0.1, 2000,
    // Numvar 500, DipRate = 0.1, Uninformative = 0.1.
    #[test]
    fn test_large_linear() {
        let config = HapCutConfig::default();
        let read_params = (2000, 10, 0.1);
        let var_params = (500, 0.1, 0.1, 0);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 7 / 10 <= largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    // Test large dataset graph.
    // HapCut with length of 10, error 0.1, 2000,
    // Numvar 500, DipRate = 0.1, Uninformative = 0.1, three loops.
    #[test]
    fn test_large_graph() {
        let config = HapCutConfig::default();
        let read_params = (2000, 10, 0.1);
        let var_params = (500, 0.1, 0.1, 3);
        let seed = 389450;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 7 / 10 <= largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    #[test]
    fn test_hap_specific_region() {
        let read_params = (100, 9, 0.1);
        let var_params = (20, 0f64, 0f64, 0);
        // There's no collapsing.
        let mut reads = gen_test_dataset(read_params, var_params, 438290);
        // Lets try to remove 3-9 in cluster 1.
        for &mut (_, ref mut read) in reads.iter_mut() {
            read.retain(|&(node, cl)| !(cl == 1 && (3..10).contains(&node)));
        }
        let copy_number: HashMap<_, _> = (0..var_params.0)
            .map(|i| match i {
                3..=9 => (i as u64, 1),
                _ => (i as u64, 2),
            })
            .collect();
        let mut config = HapCutConfig::default();
        config.copy_number = copy_number;
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
            eprintln!("{:?}", phase);
        }
        // let largest_phase = *count_phase.values().max().unwrap() as usize;
        // assert!(
        //     var_params.0 * 7 / 10 <= largest_phase,
        //     "{},{}",
        //     var_params.0,
        //     largest_phase
        // );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    #[test]
    fn test_unreliable_region() {
        let read_params = (1200, 5, 0.1);
        let var_params = (200, 0f64, 0.8f64, 0);
        let reads = gen_test_dataset(read_params, var_params, 438290);
        let config = HapCutConfig::default();
        let (phases, _consis) = hapcut(&reads, &config);
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    #[test]
    fn test_multicopy_region() {
        let read_params = (100, 10, 0.1);
        let var_params = (20, 0f64, 0f64, 0);
        // There's no collapsing.
        let mut reads = gen_test_dataset(read_params, var_params, 438290);
        // Lets try to merge 2,4,6,8,10 -> 2. It's 5 times copy.
        for &mut (_, ref mut read) in reads.iter_mut() {
            for node in read.iter_mut() {
                match node.0 {
                    2 => node.1 = node.1,
                    4 => node.1 = node.1 + 2,
                    6 => node.1 = node.1 + 4,
                    8 => node.1 = node.1 + 6,
                    _ => {}
                };
            }
        }
        let config = HapCutConfig::default();
        let (phases, _consis) = hapcut(&reads, &config);
        let mut count_phase: HashMap<_, u32> = HashMap::new();
        for phase in phases.iter() {
            *count_phase.entry(phase.phase.clone()).or_default() += 1;
        }
        let largest_phase = *count_phase.values().max().unwrap() as usize;
        assert!(
            var_params.0 * 7 / 10 <= largest_phase,
            "{},{}",
            var_params.0,
            largest_phase
        );
        let phased_blocks = split_phase_blocks(&phases, &reads, &config);
        assert_eq!(phased_blocks.len(), 2);
        let assignments: Vec<_> = reads_into_phased_block(&reads, &phased_blocks, &config)
            .into_iter()
            .map(|asn| asn.cluster)
            .collect();
        let answer = vec![vec![0; read_params.0], vec![1; read_params.0]].concat();
        let num_cl: HashSet<_> = assignments.iter().copied().collect();
        assert_eq!(num_cl.len(), 2);
        let correct = assignments
            .iter()
            .zip(answer.iter())
            .map(|(p, a)| p == a)
            .count();
        let correct = (answer.len() - correct).max(correct);
        assert!(
            90 * 2 * read_params.0 / 100 < correct,
            "{}/{}",
            correct,
            answer.len()
        );
    }
    #[test]
    fn reads_into_phased_block_test() {
        let mut phases = vec![];
        let h1: HashSet<_> = vec![(0, 0), (1, 0), (2, 0), (3, 0)].into_iter().collect();
        phases.push(h1);
        let h2: HashSet<_> = vec![(0, 1), (1, 1), (2, 1), (3, 1)].into_iter().collect();
        phases.push(h2);
        let c = HapCutConfig::default();
        assert_eq!(read_into_phased_block(&[(0, 0)], &phases, &c), Some(0));
        assert_eq!(read_into_phased_block(&[(3, 1)], &phases, &c), Some(1));
        assert_eq!(read_into_phased_block(&[(8, 1)], &phases, &c), None);
        assert_eq!(
            read_into_phased_block(&[(2, 1), (4, 0)], &phases, &c),
            Some(1)
        );
        assert_eq!(
            read_into_phased_block(&[(0, 0), (1, 1), (2, 0)], &phases, &c),
            Some(0)
        );
    }
    #[test]
    fn test_consistency() {
        let phases = vec![
            Phase::new(0, vec![-1, 1]),
            Phase::new(1, vec![-1, 1]),
            Phase::new(2, vec![-1, 1]),
            Phase::new(3, vec![-1, 1, 0]),
            Phase::new(4, vec![-1, 1, -1, 1]),
        ];
        let c = HapCutConfig::default();
        let read = vec![(0, 1), (1, 1), (2, 1)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 0);
        let read = vec![(0, 0), (1, 0), (2, 0), (4, 0), (4, 2)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 0);
        let read = vec![(0, 0), (1, 0), (2, 0), (4, 0), (4, 2)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 0);
        let read = vec![(0, 0), (1, 1), (2, 0)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 2);
        let read = vec![(0, 1), (1, 0), (2, 1)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 2);
        let read = vec![(3, 2), (3, 0), (4, 1)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 3);
        let read = vec![(3, 2), (3, 1), (4, 0), (4, 1)];
        let consistency = Phase::consistency_of_read(&phases, &read, &c);
        assert_eq!(consistency, 3);
    }
    #[test]
    fn test_weight_vector() {
        let phases = vec![
            Phase::new(0, vec![-1, 1]),
            Phase::new(1, vec![-1, 1]),
            Phase::new(2, vec![1, -1]),
            Phase::new(3, vec![-1, 1, 0]),
            Phase::new(4, vec![-1, 1, -1, 1]),
        ];
        let reads = vec![
            (0, vec![(0, 1), (1, 1), (2, 1)]),
            (1, vec![(0, 0), (1, 0), (2, 0)]),
        ];
        let weights = Phase::get_weight_vector(&phases[2], &phases, &reads);
        let answer = vec![-2f64 / 2f64, 2f64 / 2f64];
        for (w, a) in weights.0.iter().zip(answer) {
            assert!((w - a).abs() < 0.0001);
        }
        let reads = vec![
            (0, vec![(0, 1), (1, 1), (2, 1)]),
            (1, vec![(0, 0), (1, 0), (2, 0)]),
            (2, vec![(0, 1), (1, 1), (2, 0), (3, 2)]),
        ];
        let weights = Phase::get_weight_vector(&phases[2], &phases, &reads);
        let answer = vec![-1f64 + 2f64 / 3f64, 2f64 / 2f64];
        for (w, a) in weights.0.iter().zip(answer) {
            assert!((w - a).abs() < 0.0001);
        }
        let reads = vec![
            (0, vec![(0, 1), (1, 1), (2, 1)]),
            (1, vec![(0, 0), (1, 0), (2, 0)]),
            (2, vec![(0, 1), (1, 1), (2, 0), (3, 2)]),
            (3, vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (4, 2)]),
            (3, vec![(3, 0), (4, 0), (4, 3)]),
        ];
        let weights = Phase::get_weight_vector(&phases[4], &phases, &reads);
        let answer = vec![0f64, 1f64 / 5f64, 3f64 / 5f64, -1f64];
        for (w, a) in weights.0.iter().zip(answer) {
            assert!((w - a).abs() < 0.0001);
        }
    }
    #[test]
    fn test_approx() {
        let c = HapCutConfig::default();
        let target = vec![-10f64, 10f64];
        let pv = Phase::approx_parallel_vector(&target, &c);
        let answer = vec![-1, 1];
        assert_eq!(pv, answer);
        let target = vec![10f64, 10f64];
        let pv = Phase::approx_parallel_vector(&target, &c);
        let answer = vec![1, 1];
        assert_eq!(pv, answer);
        let target = vec![10f64, -10f64, 10f64];
        let pv = Phase::approx_parallel_vector(&target, &c);
        let answer = vec![1, -1, 1];
        assert_eq!(pv, answer);
        let target = vec![10f64, -0.24f64, 10f64];
        let pv = Phase::approx_parallel_vector(&target, &c);
        let answer = vec![1, 0, 1];
        assert_eq!(pv, answer);
        let target = vec![100f64, 0.15f64, 0.1f64];
        let pv = Phase::approx_parallel_vector(&target, &c);
        let answer = vec![1, 0, 0];
        assert_eq!(pv, answer);
    }
    #[test]
    fn test_phase_flip() {
        let config = HapCutConfig::default();
        let read_params = (1000, 4, 0.1);
        let var_params = (50, 0.1, 0.1, 2);
        let seed = 38940;
        let reads = gen_test_dataset(read_params, var_params, seed);
        let (mut phases, consis) = hapcut(&reads, &config);
        let changes = Phase::consistency_change_by_flip(&phases, &reads, &config);
        for (i, &change) in changes.iter().enumerate() {
            phases[i].flip();
            let consis_flip = Phase::consistency(&phases, &reads, &config);
            assert_eq!(
                consis_flip - consis,
                change,
                "{}-{}!={}",
                consis_flip,
                consis,
                change
            );
            phases[i].flip();
        }
    }
}

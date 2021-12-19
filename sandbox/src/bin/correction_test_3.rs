#[macro_use]
extern crate log;
use dirichlet_fit::Optimizer;
use rand::{Rng, SeedableRng};
use rand_distr::{self, Distribution};
use rand_xoshiro::Xoshiro256PlusPlus;
// If the posterior prob is less than this value, it would be removed...
const WEIGHT_FILTER: f64 = 0.01;
// (id, [(unit, posterior distr)]
type Haplotype = (Vec<(usize, usize)>, usize, Vec<(usize, usize)>);
type ClusterNum = (Vec<usize>, usize, Vec<usize>);
type IsError = (Vec<bool>, bool, Vec<bool>);
type UnitArray = Vec<(usize, Vec<f64>)>;
type SimRead = (UnitArray, Vec<f64>, UnitArray);
fn sim_read<R: Rng>(
    (hap_up, hap, hap_down): &Haplotype,
    (up_cl_num, cl_num, down_cl_num): &ClusterNum,
    (up_err, err, down_err): &IsError,
    rng: &mut R,
) -> SimRead {
    let upstream: Vec<_> = hap_up
        .iter()
        .zip(up_cl_num.iter())
        .zip(up_err.iter())
        .map(|((&(unit, cl), &dim), &err)| (unit, gen_post_dist(cl, dim, err, rng)))
        .collect();
    let center = gen_post_dist(*hap, *cl_num, *err, rng);
    let downstream: Vec<_> = hap_down
        .iter()
        .zip(down_cl_num.iter())
        .zip(down_err.iter())
        .map(|((&(unit, cl), &dim), &err)| (unit, gen_post_dist(cl, dim, err, rng)))
        .collect();
    (upstream, center, downstream)
}

fn gen_post_dist<R: Rng>(cl: usize, dim: usize, err: bool, rng: &mut R) -> Vec<f64> {
    let mut param = vec![0.5f64; dim];
    if !err {
        param[cl] = 20.0;
    }
    let dirichlet = rand_distr::Dirichlet::new(&param).unwrap();
    dirichlet.sample(rng).into_iter().map(|x| x.ln()).collect()
}

fn gen_reads<R: Rng>(
    hap1: &Haplotype,
    hap2: &Haplotype,
    clnum: &ClusterNum,
    errors: &IsError,
    rng: &mut R,
    num: usize,
) -> Vec<Context> {
    (0..2 * num)
        .map(|i| match i / num {
            0 => (i, sim_read(hap1, clnum, errors, rng)),
            1 => (i, sim_read(hap2, clnum, errors, rng)),
            _ => panic!(),
        })
        .map(|(id, (upstream, center, downstream))| Context {
            id: id as u64,
            index: 0,
            unit: 0,
            center,
            upstream,
            downstream,
        })
        .collect()
}

fn main() {
    env_logger::init();
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(4234);
    // let hap1s = {
    //     let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
    //     let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
    //     (hap1_up, 0, hap1_down)
    // };
    // let hap2s = {
    //     let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
    //     let hap2_down = vec![(0, 1), (1, 1), (2, 1)];

    //     (hap2_up, 0, hap2_down)
    // };
    // let cl_num = { (vec![2, 2, 2], 2, vec![2, 2, 2]) };
    // let is_uninformative = { (vec![false, false, false], false, vec![false, false, false]) };
    let hap1s = (vec![], 0, vec![]);
    let hap2s = (vec![], 1, vec![]);
    let cl_num = { (vec![], 2, vec![]) };
    let is_uninformative = (vec![], false, vec![]);
    let num = 10;
    let reads = gen_reads(&hap1s, &hap2s, &cl_num, &is_uninformative, &mut rng, num);
    for (i, read) in reads.iter().enumerate() {
        trace!("{}\t{}", i, vec2str(&read.center));
    }
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(4);
    let (result, _) = clustering(&reads, 2, &mut rng);
    for (id, _, post) in result {
        eprintln!("{}\t{:?}", id, post);
    }
}

#[derive(Debug, Clone)]
struct Context {
    // The ID of the original encoded read.
    id: u64,
    // The original index of this context.
    index: usize,
    #[allow(dead_code)]
    unit: u64,
    // This is log-scaled.
    center: Vec<f64>,
    // Only forward/backward information is preversed, the
    // ordering is arbitrary.
    // LOG posterior probabilities.
    upstream: Vec<(usize, Vec<f64>)>,
    downstream: Vec<(usize, Vec<f64>)>,
}

fn clustering<R: Rng>(
    contexts: &[Context],
    k: usize,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64) {
    // ID of this trial.
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    let mut weights = if k == 1 {
        vec![vec![1f64]; contexts.len()]
    } else {
        let dir = rand_distr::Dirichlet::new(&vec![0.5f64; k]).unwrap();
        contexts.iter().map(|_| dir.sample(rng)).collect()
    };
    for (weights, ctx) in weights.iter().zip(contexts.iter()) {
        trace!("{}\t{}", vec2str(&ctx.center), vec2str(weights));
    }
    let mut model = HMMixtureModel::new(contexts, &weights, k);
    trace!("CORRECT\tModel\t{}\t{}\n{}", -1, id, model);
    let mut lk: f64 = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
    trace!("CORRECT\tLikelihood\t{}\t{}", id, lk);
    for t in 0..20 {
        trace!("CORRECT\tModel\t{}\t{}\n{}", t, id, model);
        model.update(&mut weights, contexts);
        let next_lk = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
        trace!("CORRECT\tLikelihood\t{}\t{}", id, next_lk);
        if (next_lk - lk) < 0.001 {
            break;
        } else {
            lk = next_lk;
        }
    }
    trace!("CORRECT\tModel\t{}\n{}", id, model);
    let predictions: Vec<_> = contexts
        .iter()
        .zip(weights.into_iter())
        .map(|(ctx, weight)| (ctx.id, ctx.index, weight))
        .collect();
    trace!("CORRECT\tFinal\t{}\t{}\t{:.4}", id, k, lk);
    (predictions, lk)
}

struct HMMixtureModel {
    k: usize,
    // Number of upstream units.
    upstream_units: usize,
    // Number of dowstream units.
    downstream_units: usize,
    fractions: Vec<f64>,
    models: Vec<HMModel>,
}

impl std::fmt::Display for HMMixtureModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}\t{}\t{}",
            self.k, self.upstream_units, self.downstream_units
        )?;
        for (i, (fr, m)) in self.fractions.iter().zip(self.models.iter()).enumerate() {
            if i != self.k - 1 {
                writeln!(f, "{}\t{:.3}\t{}", i, fr, m)?;
            } else {
                write!(f, "{}\t{:.3}\t{}", i, fr, m)?;
            }
        }
        Ok(())
    }
}

impl HMMixtureModel {
    fn new(contexts: &[Context], weights: &[Vec<f64>], cluster_num: usize) -> Self {
        let fractions = sum_and_normalize(weights);
        let up_units = contexts
            .iter()
            .filter_map(|ctx| ctx.upstream.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap_or(0);
        let down_units = contexts
            .iter()
            .filter_map(|ctx| ctx.downstream.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap_or(0);
        let models: Vec<_> = (0..cluster_num)
            .map(|cl| HMModel::new(contexts, weights, cl, up_units, down_units))
            .collect();
        Self {
            k: cluster_num,
            upstream_units: up_units,
            downstream_units: down_units,
            fractions,
            models,
        }
    }
    fn get_likelihood(&self, context: &Context) -> f64 {
        let lks: Vec<_> = self
            .fractions
            .iter()
            .zip(self.models.iter())
            .map(|(f, m)| f.ln() + m.lk(context))
            .collect();
        logsumexp(&lks)
    }
    fn update(&mut self, weights: &mut [Vec<f64>], contexts: &[Context]) {
        // E-step
        // Posterior prob and alignment expectation.
        let (posterior_prob, align_expt) = self.e_step(contexts);
        for (i, post) in posterior_prob.iter().enumerate() {
            trace!("{}\t[{}]", i, vec2str(post));
        }
        // M-step
        self.fractions = sum_and_normalize(&posterior_prob);
        for (cl, model) in self.models.iter_mut().enumerate() {
            // We should filter out elements with very small weight.
            let (mut ws, mut ctces, mut alns) = (vec![], vec![], vec![]);
            for ((weight, ctx), align) in weights
                .iter()
                .zip(contexts.iter())
                .zip(align_expt[cl].iter())
                .filter(|&((ws, _), _)| 0.01 < ws[cl])
            {
                ws.push(weight[cl]);
                ctces.push(ctx);
                alns.push(align);
            }
            assert_eq!(ws.len(), ctces.len());
            assert_eq!(ws.len(), alns.len());
            model.update(&ws, &ctces, &alns);
        }
        for (ws, posterior) in weights.iter_mut().zip(posterior_prob) {
            *ws = posterior;
        }
    }
    // 2nd return value: the expectation of aligning information of each cluster k.
    fn e_step(&self, contexts: &[Context]) -> (Vec<Vec<f64>>, Vec<Vec<AlignInfo>>) {
        let alignments: Vec<Vec<_>> = self
            .models
            .iter()
            .map(|m| contexts.iter().map(|ctx| m.align(ctx)).collect())
            .collect();
        let posterior: Vec<_> = (0..contexts.len())
            .map(|idx| {
                let mut lks: Vec<_> = self
                    .fractions
                    .iter()
                    .zip(alignments.iter())
                    .map(|(f, alns)| f.ln() + alns[idx].lk)
                    .collect();
                let total = logsumexp(&lks);
                lks.iter_mut().for_each(|x| *x = (*x - total).exp());
                lks
            })
            .collect();
        (posterior, alignments)
    }
}

fn vec2str(xs: &[f64]) -> String {
    let xs: Vec<_> = xs.iter().map(|x| format!("{:.2}", x)).collect();
    xs.join(",")
}

struct AlignInfo {
    // Likelihood. Sum of center, upstream and downstream.
    lk: f64,
    // For each j-th position of the reference and the i-th position of the query,
    // upstream[j][i] return the probability that the j-th pos of the refr emit the i-th pos of the query.
    // Note that the last element of each row(upstream[j].last()) would be the prob of deletion.
    upstream: Vec<Vec<f64>>,
    downstream: Vec<Vec<f64>>,
}

// Sum xs and normalize the result so that it summing up to 1.
fn sum_and_normalize(xss: &[Vec<f64>]) -> Vec<f64> {
    let mut sumed = vec![0f64; xss[0].len()];
    let mut total = 0f64;
    for xs in xss {
        sumed
            .iter_mut()
            .zip(xs.iter())
            .for_each(|(acc, x)| *acc += x);
        total += xs.iter().sum::<f64>();
    }
    sumed.iter_mut().for_each(|x| *x /= total);
    sumed
}

impl std::fmt::Display for HMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for up in self.upstream.iter() {
            writeln!(f, "{}", up)?;
        }
        writeln!(f, "{}", self.center)?;
        for (i, down) in self.downstream.iter().enumerate() {
            if i != self.downstream_len - 1 {
                writeln!(f, "{}", down)?;
            } else {
                write!(f, "{}", down)?;
            }
        }
        Ok(())
    }
}

struct HMModel {
    upstream_len: usize,
    downstream_len: usize,
    center: Dirichlet,
    // the f64 parts summing up to 1 on each position.
    upstream: Vec<DirichletMixture>,
    // the f64 parts summing up to 1 on each position.
    downstream: Vec<DirichletMixture>,
}

const NEG_LARGE: f64 = -10000000000000f64;
type DP = Vec<Vec<f64>>;
impl HMModel {
    fn new(
        contexts: &[Context],
        weights: &[Vec<f64>],
        cluster: usize,
        up: usize,
        down: usize,
    ) -> Self {
        let upstream_len: usize = contexts
            .iter()
            .map(|ctx| ctx.upstream.len())
            .max()
            .unwrap_or(0);
        let downstream_len: usize = contexts
            .iter()
            .map(|ctx| ctx.downstream.len())
            .max()
            .unwrap_or(0);
        let weights: Vec<_> = weights.iter().map(|ws| ws[cluster]).collect();
        let center: Vec<_> = contexts.iter().map(|ctx| ctx.center.as_slice()).collect();
        let center = Dirichlet::new(&center, &weights);
        let mut upstream_slots = vec![vec![]; upstream_len];
        for (ctx, &w) in contexts.iter().zip(weights.iter()) {
            for (&(unit, ref prob), slot) in ctx.upstream.iter().zip(upstream_slots.iter_mut()) {
                slot.push((unit, w, prob.as_slice()));
            }
        }
        let upstream: Vec<_> = upstream_slots
            .iter()
            .map(|bucket| DirichletMixture::new(bucket, up))
            .collect();
        let mut downstream_slots = vec![vec![]; downstream_len];
        for (ctx, &w) in contexts.iter().zip(weights.iter()) {
            for (&(unit, ref prob), slot) in ctx.downstream.iter().zip(downstream_slots.iter_mut())
            {
                slot.push((unit, w, prob.as_slice()));
            }
        }
        let downstream: Vec<_> = downstream_slots
            .iter()
            .map(|bucket| DirichletMixture::new(bucket, down))
            .collect();
        Self {
            upstream_len,
            downstream_len,
            center,
            upstream,
            downstream,
        }
    }
    fn lk(&self, ctx: &Context) -> f64 {
        let up_forward = Self::forward(&self.upstream, &ctx.upstream);
        let up_lk = up_forward[ctx.upstream.len()][self.upstream_len];
        let down_forward = Self::forward(&self.downstream, &ctx.downstream);
        let down_lk = down_forward[ctx.downstream.len()][self.downstream_len];
        let center = self.center.lk(&ctx.center);
        up_lk + center + down_lk
    }
    fn align(&self, ctx: &Context) -> AlignInfo {
        let up_forward = Self::forward(&self.upstream, &ctx.upstream);
        let up_backward = Self::backward(&self.upstream, &ctx.upstream);
        let down_forward = Self::forward(&self.downstream, &ctx.downstream);
        let down_backward = Self::backward(&self.downstream, &ctx.downstream);
        let up_lk = up_forward[ctx.upstream.len()][self.upstream_len];
        let down_lk = down_forward[ctx.downstream.len()][self.downstream_len];
        let center = self.center.lk(&ctx.center);
        let upstream = Self::to_align(&self.upstream, &ctx.upstream, &up_forward, &up_backward);
        let downstream = Self::to_align(
            &self.downstream,
            &ctx.downstream,
            &down_forward,
            &down_backward,
        );
        AlignInfo {
            lk: up_lk + center + down_lk,
            upstream,
            downstream,
        }
    }
    fn forward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)]) -> DP {
        // Init
        let mut dp = vec![vec![0f64; query.len() + 1]; refr.len() + 1];
        for i in 1..query.len() + 1 {
            dp[i][0] = NEG_LARGE;
        }
        for j in 1..refr.len() + 1 {
            dp[0][j] = dp[0][j - 1] + refr[j - 1].del();
        }
        for (i, obs) in query.iter().enumerate().map(|(i, x)| (i + 1, x)) {
            for (j, dir) in refr.iter().enumerate().map(|(j, x)| (j + 1, x)) {
                let del_trans = dp[i][j - 1] + dir.del();
                let match_trans = dp[i - 1][j - 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans);
            }
        }
        dp
    }
    fn backward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)]) -> DP {
        // Init
        let mut dp = vec![vec![0f64; query.len() + 1]; refr.len() + 1];
        for i in 1..query.len() + 1 {
            dp[i][refr.len()] = NEG_LARGE;
        }
        for j in 1..refr.len() + 1 {
            dp[query.len()][j] = 0f64;
        }
        for (i, obs) in query.iter().enumerate().rev() {
            for (j, dir) in refr.iter().enumerate().rev() {
                let del_trans = dp[i][j + 1] + dir.del();
                let match_trans = dp[i + 1][j + 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans);
            }
        }
        dp
    }
    fn to_align(
        refr: &[DirichletMixture],
        query: &[(usize, Vec<f64>)],
        forward: &DP,
        backward: &DP,
    ) -> Vec<Vec<f64>> {
        refr.iter()
            .enumerate()
            .map(|(j, dir)| {
                let mut lks: Vec<_> = query
                    .iter()
                    .enumerate()
                    .map(|(i, obs)| forward[i][j] + dir.mat(obs) + backward[i + 1][j + 1])
                    .collect();
                // Push deletion
                {
                    let del_probs: Vec<_> = (0..query.len())
                        .map(|i| forward[i][j] + dir.del() + backward[i][j + 1])
                        .collect();
                    lks.push(logsumexp(&del_probs));
                }
                let total = logsumexp(&lks);
                lks.iter_mut().for_each(|x| *x = (*x - total).exp());
                lks
            })
            .collect()
    }
    fn update(&mut self, weights: &[f64], contexts: &[&Context], align_expt: &[&AlignInfo]) {
        let center: Vec<_> = contexts.iter().map(|ctx| ctx.center.as_slice()).collect();
        self.center.update(&center, weights);
        // Update upstream
        let up_align_expectation: Vec<_> = align_expt
            .iter()
            .map(|aln| aln.upstream.as_slice())
            .collect();
        let up_contexts: Vec<_> = contexts.iter().map(|ctx| ctx.upstream.as_slice()).collect();
        Self::update_oneside(
            &mut self.upstream,
            weights,
            &up_contexts,
            &up_align_expectation,
        );
        let down_align_expectation: Vec<&[Vec<f64>]> = align_expt
            .iter()
            .map(|aln| aln.downstream.as_slice())
            .collect();
        let down_contexts: Vec<_> = contexts
            .iter()
            .map(|ctx| ctx.downstream.as_slice())
            .collect();
        Self::update_oneside(
            &mut self.downstream,
            weights,
            &down_contexts,
            &down_align_expectation,
        );
    }
    fn update_oneside(
        parameters: &mut [DirichletMixture],
        weights: &[f64],
        contexts: &[&[(usize, Vec<f64>)]],
        alignments: &[&[Vec<f64>]],
    ) {
        for (j, dir) in parameters.iter_mut().enumerate() {
            let alns = alignments.iter().map(|aln| aln[j].iter());
            let total_ws: f64 = alns
                .clone()
                .zip(weights.iter())
                .map(|(alns, w)| w * alns.sum::<f64>())
                .sum();
            // It never panics.
            let del_ws: f64 = alns
                .clone()
                .zip(weights.iter())
                .map(|(alns, w)| w * alns.last().unwrap())
                .sum();
            dir.del_prob = del_ws / total_ws;
            // First re-estimate parameters on categorical distributions.
            let mut fraction = vec![0f64; dir.dirichlets.len()];
            for ((w, ctx), alns) in weights.iter().zip(contexts.iter()).zip(alns.clone()) {
                for (&(unit, _), weight) in ctx.iter().zip(alns) {
                    fraction[unit] += w * weight;
                }
            }
            let total: f64 = fraction.iter().sum();
            dir.dirichlets
                .iter_mut()
                .zip(fraction)
                .for_each(|(x, f)| x.0 = f / total);
            // Then, re-estimate parameters on dirichlet distirbution on high-frequent unit.
            for (target, (_, dir)) in dir
                .dirichlets
                .iter_mut()
                .enumerate()
                .filter(|(_, (fr, dir))| WEIGHT_FILTER < *fr && 1 < dir.dim)
            {
                let (weights, dataset): (Vec<_>, Vec<_>) = weights
                    .iter()
                    .zip(contexts.iter())
                    .zip(alns.clone())
                    .map(|((&w, ctx), alns)| {
                        let (ws, obss): (Vec<f64>, Vec<_>) = ctx
                            .iter()
                            .zip(alns)
                            .filter_map(|(&(unit, ref prob), &weight)| {
                                (unit == target).then(|| (weight, prob.as_slice()))
                            })
                            .unzip();
                        ((w, ws), obss)
                    })
                    .unzip();
                assert!(dir.dim > 1);
                let mut opt = dirichlet_fit::AdamOptimizer::new(dir.dim);
                dir.param =
                    dirichlet_fit::fit_multiple_with(&dataset, &weights, &mut opt, &dir.param);
            }
        }
    }
}

struct DirichletMixture {
    del_prob: f64,
    dirichlets: Vec<(f64, Dirichlet)>,
}

impl std::fmt::Display for DirichletMixture {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dump: Vec<_> = self
            .dirichlets
            .iter()
            .enumerate()
            .filter(|&(_, &(fr, _))| fr > 0.1)
            .map(|(unit, (fr, dir))| format!("{}:{:.3}:{}", unit, fr, dir))
            .collect();
        write!(f, "{}\t{}", self.del_prob, dump.join("\t"))
    }
}

impl DirichletMixture {
    // Return log del_prob
    fn del(&self) -> f64 {
        self.del_prob.ln()
    }
    // Return log (1-e)Dir(prob|param).
    fn mat(&self, &(unit, ref prob): &(usize, Vec<f64>)) -> f64 {
        let &(frac, ref dir) = &self.dirichlets[unit];
        (1f64 - self.del_prob).ln() + frac.ln() + dir.lk(prob)
    }
    // Create new instance from the given data.
    fn new(observed: &[(usize, f64, &[f64])], len: usize) -> Self {
        let del_prob = 0.1;
        let mut post_probs = vec![vec![]; len];
        let mut prob_weights = vec![vec![]; len];
        // TODO: Maybe 1 would be better?
        let mut fractions = vec![0f64; len];
        for &(unit, w, prob) in observed {
            fractions[unit] += w;
            post_probs[unit].push(prob);
            prob_weights[unit].push(w);
        }
        let sum: f64 = fractions.iter().sum();
        let dirichlets: Vec<_> = post_probs
            .iter()
            .zip(prob_weights.iter())
            .zip(fractions.iter())
            .map(|((probs, ws), fr)| {
                assert_eq!(probs.len(), ws.len());
                (fr / sum, Dirichlet::new(probs, ws))
            })
            .collect();
        Self {
            del_prob,
            dirichlets,
        }
    }
}

struct Dirichlet {
    dim: usize,
    param: Vec<f64>,
}

impl Dirichlet {
    fn lk(&self, log_prob: &[f64]) -> f64 {
        assert_eq!(self.dim, log_prob.len());
        match self.dim {
            1 => 0f64,
            _ => dirichlet_fit::dirichlet_log(log_prob, &self.param),
        }
    }
    fn update<T: std::borrow::Borrow<[f64]>>(&mut self, obs: &[T], weights: &[f64]) {
        if 1 < self.dim {
            let (data, weights) = (&[obs], &[(1f64, weights)]);
            let mut optimizer = dirichlet_fit::AdamOptimizer::new(self.dim);
            self.param =
                dirichlet_fit::fit_multiple_with(data, weights, &mut optimizer, &self.param);
        }
    }
    fn new<T: std::borrow::Borrow<[f64]>>(obs: &[T], weights: &[f64]) -> Self {
        let dim = obs[0].borrow().len();
        let param = match dim {
            1 => vec![1f64],
            _ => dirichlet_fit::fit_multiple(&[obs], &[(1f64, weights)]),
        };
        Self { dim, param }
    }
}

impl std::fmt::Display for Dirichlet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, p) in self.param.iter().enumerate() {
            if i != self.dim - 1 {
                write!(f, "{:.2},", p)?;
            } else {
                write!(f, "{:.2}", p)?;
            }
        }
        write!(f, "]")
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        0.
    } else if xs.len() == 1 {
        xs[0]
    } else {
        let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
        assert!(sum >= 0., "{:?}->{}", xs, sum);
        max + sum
    }
}
// Log(exp(x) + exp(y)) = x + Log(1+exp(y-x)).
fn logsumexp2(x: f64, y: f64) -> f64 {
    let (x, y) = (x.max(y), x.min(y));
    x + (1f64 + (y - x).exp()).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lse() {
        let x = 10.3f64;
        let y = 3.1f64;
        let test = logsumexp2(x, y);
        let answer = (x.exp() + y.exp()).ln();
        assert!((test - answer).abs() < 0.0001);
    }
}

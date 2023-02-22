use rand::{prelude::SliceRandom, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
pub fn estimate_minimum_gain(hmm: &PairHiddenMarkovModelOnStrands) -> f64 {
    const SEED: u64 = 23908;
    const SAMPLE_NUM: usize = 1000;
    const SEQ_NUM: usize = 500;
    const LEN: usize = 100;
    const BAND: usize = 25;
    const MIN_REQ: f64 = 1f64;
    let mut medians: Vec<_> = (0..SAMPLE_NUM)
        .into_par_iter()
        .map(|seed| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(SEED + seed as u64);
            let hap1 = kiley::gen_seq::generate_seq(&mut rng, LEN);
            let hap2 = kiley::gen_seq::introduce_errors(&hap1, &mut rng, 0, 1, 0);
            use kiley::gen_seq::Generate;
            let mut lks: Vec<_> = (0..SEQ_NUM)
                .map(|t| {
                    let hmm = match t % 2 == 0 {
                        true => hmm.forward(),
                        false => hmm.reverse(),
                    };
                    let read = hmm.gen(&hap1, &mut rng);
                    let lk_base = hmm.likelihood_antidiagonal_bootstrap(&hap1, &read, BAND);
                    let lk_diff = hmm.likelihood_antidiagonal_bootstrap(&hap2, &read, BAND);
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

#[derive(Debug, Clone, Copy)]
pub struct GainProfile {
    // Median gain of this variant.
    gain: f64,
    // Probability to see this variant from the null model.
    prob: f64,
}

impl std::fmt::Display for GainProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.0},{:.2}", self.gain, self.prob)
    }
}

#[derive(Debug, Clone)]
pub struct Gains {
    max_homopolymer_len: usize,
    subst: Vec<GainProfile>,
    deletions: Vec<GainProfile>,
    insertions: Vec<GainProfile>,
}

impl std::fmt::Display for Gains {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let subst: Vec<_> = self.subst.iter().map(|x| format!("{x}")).collect();
        let deletions: Vec<_> = self.deletions.iter().map(|x| format!("{x}")).collect();
        let insertions: Vec<_> = self.insertions.iter().map(|x| format!("{x}")).collect();
        write!(
            f,
            "{}\n{}\n{}",
            subst.join("\t"),
            deletions.join("\t"),
            insertions.join("\t")
        )
    }
}

impl Gains {
    pub fn expected(&self, homop_len: usize, diff_type: DiffType) -> f64 {
        assert!(0 < homop_len);
        let homop_len = homop_len.min(self.max_homopolymer_len);
        match diff_type {
            DiffType::Subst => self.subst[homop_len - 1].gain,
            DiffType::Del => self.deletions[homop_len - 1].gain,
            DiffType::Ins => self.insertions[homop_len - 1].gain,
        }
    }
    pub fn pvalues(&self, total: usize) -> Pvalues {
        let substs: Vec<_> = self
            .subst
            .iter()
            .map(|gp| pvalues(gp.prob, total))
            .collect();
        let deletions: Vec<_> = self
            .deletions
            .iter()
            .map(|gp| pvalues(gp.prob, total))
            .collect();
        let insertions: Vec<_> = self
            .insertions
            .iter()
            .map(|gp| pvalues(gp.prob, total))
            .collect();
        Pvalues {
            max_homopolymer_len: self.max_homopolymer_len,
            total,
            substs,
            deletions,
            insertions,
        }
    }
}

// i -> Prob(i <= X | n, prob)
fn pvalues(prob: f64, n: usize) -> Vec<f64> {
    let (ln, in_ln) = (prob.ln(), (1f64 - prob).ln());
    let mut log_probs = vec![in_ln * n as f64];
    for k in 0..n {
        let prev = *log_probs.last().unwrap();
        let offset = ln + ((n - k) as f64).ln() - in_ln - ((k + 1) as f64).ln();
        log_probs.push(prev + offset);
    }
    let mut cum_sum = log_probs;
    for k in (0..n).rev() {
        cum_sum[k] = logsumexp(cum_sum[k + 1], cum_sum[k]);
    }
    cum_sum.iter_mut().for_each(|x| *x = x.exp());
    cum_sum
}

fn logsumexp(x: f64, y: f64) -> f64 {
    if y < x {
        x + (1f64 + (y - x).exp()).ln()
    } else {
        y + (1f64 + (x - y).exp()).ln()
    }
}

#[derive(Debug, Clone)]
pub struct Pvalues {
    max_homopolymer_len: usize,
    total: usize,
    substs: Vec<Vec<f64>>,
    deletions: Vec<Vec<f64>>,
    insertions: Vec<Vec<f64>>,
}

impl Pvalues {
    pub fn pvalue(&self, homop_len: usize, diff_type: DiffType, count: usize) -> f64 {
        assert!(0 < homop_len);
        assert!(count <= self.total);
        let homop_len = homop_len.min(self.max_homopolymer_len);
        match diff_type {
            DiffType::Subst => self.substs[homop_len - 1][count],
            DiffType::Del => self.deletions[homop_len - 1][count],
            DiffType::Ins => self.insertions[homop_len - 1][count],
        }
    }
}

use kiley::hmm::PairHiddenMarkovModelOnStrands;
pub fn estimate_gain(
    hmm: &PairHiddenMarkovModelOnStrands,
    seed: u64,
    seq_len: usize,
    band: usize,
    homop_len: usize,
) -> Gains {
    let subst: Vec<_> = (1..=homop_len)
        .map(|len| gain_of(hmm, seed, seq_len, band, len, DiffType::Subst))
        .collect();
    let deletions: Vec<_> = (1..=homop_len)
        .map(|len| gain_of(hmm, seed, seq_len, band, len, DiffType::Del))
        .collect();
    let insertions: Vec<_> = (1..=homop_len)
        .map(|len| gain_of(hmm, seed, seq_len, band, len, DiffType::Ins))
        .collect();
    Gains {
        max_homopolymer_len: homop_len,
        subst,
        deletions,
        insertions,
    }
}

const SEED: u64 = 309423;
const SEQ_LEN: usize = 100;
const BAND: usize = 10;
const HOMOP_LEN: usize = 3;
pub fn estimate_gain_default(hmm: &PairHiddenMarkovModelOnStrands) -> Gains {
    crate::likelihood_gains::estimate_gain(hmm, SEED, SEQ_LEN, BAND, HOMOP_LEN)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiffType {
    Subst,
    Del,
    Ins,
}

impl std::fmt::Display for DiffType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DiffType::Subst => write!(f, "S"),
            DiffType::Del => write!(f, "D"),
            DiffType::Ins => write!(f, "I"),
        }
    }
}

const BASES: &[u8] = b"ACGT";
use rand::Rng;
fn sample_triple<R: Rng>(rng: &mut R) -> (u8, u8, u8) {
    let homop = *BASES.choose(rng).unwrap();
    let right = *BASES
        .choose_weighted(rng, |&b| if b != homop { 1f64 } else { 0.0 })
        .unwrap();
    let left = *BASES
        .choose_weighted(rng, |&b| if b != homop && b != right { 1f64 } else { 0.0 })
        .unwrap();
    (right, homop, left)
}

fn gen_diff_haplotypes<R: Rng>(rng: &mut R, len: usize, diff_type: DiffType) -> (Vec<u8>, Vec<u8>) {
    let (right, center, left) = sample_triple(rng);
    let center1 = vec![center; len];
    let center2 = {
        let mut center2 = center1.clone();
        match diff_type {
            DiffType::Subst => {
                let diff = *BASES
                    .choose_weighted(rng, |&b| if b != center { 1f64 } else { 0.0 })
                    .unwrap();
                center2[0] = diff;
            }
            DiffType::Del => {
                center2.remove(0);
            }
            DiffType::Ins => {
                let diff = *BASES
                    .choose_weighted(rng, |&b| if b != center { 1f64 } else { 0.0 })
                    .unwrap();
                center2.insert(1, diff);
            }
        }
        center2
    };
    let hap1 = vec![vec![right], center1, vec![left]].concat();
    let hap2 = vec![vec![right], center2, vec![left]].concat();
    (hap1, hap2)
}

fn gain_of(
    hmm: &PairHiddenMarkovModelOnStrands,
    seed: u64,
    seq_len: usize,
    band: usize,
    len: usize,
    diff_type: DiffType,
) -> GainProfile {
    const SAMPLE_NUM: usize = 100;
    const GAIN_POS: usize = SAMPLE_NUM / 10;
    const PROB_POS: usize = SAMPLE_NUM * 2 / 3;
    const SEQ_NUM: usize = 50;
    use kiley::gen_seq::Generate;
    let (mut medians, mut probs): (Vec<_>, Vec<_>) = (0..SAMPLE_NUM)
        .into_par_iter()
        .map(|i| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(i as u64 + seed);
            let seg1 = kiley::gen_seq::generate_seq(&mut rng, seq_len / 2);
            let seg2 = kiley::gen_seq::generate_seq(&mut rng, seq_len / 2);
            let (hap1, hap2) = gen_diff_haplotypes(&mut rng, len, diff_type);
            let template = vec![seg1.clone(), hap1, seg2.clone()].concat();
            let diff = vec![seg1, hap2, seg2].concat();
            let mut lk_diff: Vec<_> = (0..SEQ_NUM)
                .map(|t| {
                    let hmm = match t % 2 == 0 {
                        true => hmm.forward(),
                        false => hmm.reverse(),
                    };
                    let read = hmm.gen(&diff, &mut rng);
                    let lk_base = hmm.likelihood_antidiagonal_bootstrap(&template, &read, band);
                    let lk_diff = hmm.likelihood_antidiagonal_bootstrap(&diff, &read, band);
                    lk_diff - lk_base
                })
                .collect();
            let (_, &mut expected_gain, _) =
                lk_diff.select_nth_unstable_by(SEQ_NUM / 2, |x, y| x.partial_cmp(y).unwrap());
            let min_gain = match diff_type {
                DiffType::Subst => expected_gain / 10f64,
                DiffType::Del => 0.0001,
                DiffType::Ins => 0.0001,
            };
            let null_prob = (0..SEQ_NUM)
                .filter(|t| {
                    let hmm = match t % 2 == 0 {
                        true => hmm.forward(),
                        false => hmm.reverse(),
                    };
                    let read = hmm.gen(&template, &mut rng);
                    let lk_base = hmm.likelihood_antidiagonal_bootstrap(&template, &read, band);
                    let lk_diff = hmm.likelihood_antidiagonal_bootstrap(&diff, &read, band);
                    lk_base + min_gain < lk_diff
                })
                .count();
            (expected_gain, null_prob as f64 / SEQ_NUM as f64)
        })
        .collect();
    let (_, &mut gain, _) =
        medians.select_nth_unstable_by(GAIN_POS, |x, y| x.partial_cmp(y).unwrap());
    let (_, &mut prob, _) =
        probs.select_nth_unstable_by(PROB_POS, |x, y| x.partial_cmp(y).unwrap());
    let prob = prob.max(0.000000001);
    GainProfile { gain, prob }
}

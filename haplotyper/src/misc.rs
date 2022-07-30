//! Useful functions.

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

/// Input: observation of each occurence,
/// return Cramer's V statistics.
pub fn cramers_v(labels: &[(u32, u32)], (cl1, cl2): (usize, usize)) -> f64 {
    let mut first_occs = vec![0; cl1];
    let mut second_occs = vec![0; cl2];
    let mut occs = vec![vec![0; cl2]; cl1];
    for &(f, s) in labels.iter() {
        first_occs[f as usize] += 1;
        second_occs[s as usize] += 1;
        occs[f as usize][s as usize] += 1;
    }
    let expected: Vec<Vec<_>> = first_occs
        .iter()
        .map(|focc| {
            second_occs
                .iter()
                .map(|socc| (focc * socc) as f64 / labels.len() as f64)
                .collect()
        })
        .collect();
    let chi_sq: f64 = expected
        .iter()
        .zip(occs.iter())
        .map(|(exp, occ)| -> f64 {
            exp.iter()
                .zip(occ.iter())
                .filter(|&(&e, _)| 0f64 < e)
                .map(|(e, &o)| (e - o as f64).powi(2) / e)
                .sum()
        })
        .sum();
    let denom = labels.len() * (cl1.min(cl2) - 1);
    assert!(0 < denom, "{:?}", labels);
    (chi_sq / denom as f64).sqrt()
}

pub fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

pub fn logsumexp_str<I: Iterator<Item = f64>>(xs: I) -> f64 {
    let (mut max, mut accum, mut count) = (std::f64::NEG_INFINITY, 0f64, 0);
    for x in xs {
        count += 1;
        if x <= max {
            accum += (x - max).exp();
        } else {
            accum = (max - x).exp() * accum + 1f64;
            max = x;
        }
    }
    match count {
        1 => max,
        _ => accum.ln() + max,
    }
}

const EDLIB2KILEY: [kiley::Op; 4] = [
    kiley::Op::Match,
    kiley::Op::Ins,
    kiley::Op::Del,
    kiley::Op::Mismatch,
];
pub fn edlib_to_kiley(edlib: &[u8]) -> Vec<kiley::Op> {
    edlib.iter().map(|&op| EDLIB2KILEY[op as usize]).collect()
}

pub fn ops_to_kiley(ops: &definitions::Ops) -> Vec<kiley::Op> {
    use definitions::Op;
    ops.iter()
        .flat_map(|op| match op {
            Op::Match(l) => std::iter::repeat(kiley::Op::Match).take(*l),
            Op::Del(l) => std::iter::repeat(kiley::Op::Del).take(*l),
            Op::Ins(l) => std::iter::repeat(kiley::Op::Ins).take(*l),
        })
        .collect()
}

use rand::Rng;
const UPDATE_THR: f64 = 0.00000001;
/// kmenas++. Returns the residual and the assginemtns.
pub fn kmeans<R: Rng, D: std::borrow::Borrow<[f64]>>(
    data: &[D],
    k: usize,
    rng: &mut R,
) -> (f64, Vec<usize>) {
    assert!(1 <= k);
    let dim = data[0].borrow().len();
    assert!(0 < dim);
    let mut assignments = match rng.gen_bool(0.5) {
        true => (0..data.len()).map(|_| rng.gen_range(0..k)).collect(),
        false => suggest_first(data, k, rng),
    };
    let (mut centers, mut counts) = (vec![vec![0f64; dim]; k], vec![0; k]);
    let mut dist = get_dist(data, &centers, &assignments);
    loop {
        // Update
        update_centers(data, &mut centers, &mut counts, &assignments);
        update_assignments(data, &centers, &mut assignments);
        let new_dist = get_dist(data, &centers, &assignments);
        assert!(new_dist <= dist);
        assert!(0f64 <= dist - new_dist);
        if dist - new_dist < UPDATE_THR {
            break;
        } else {
            dist = new_dist;
        }
    }
    (dist, assignments)
}

fn update_assignments<D: std::borrow::Borrow<[f64]>, C: std::borrow::Borrow<[f64]>>(
    data: &[D],
    centers: &[C],
    assignments: &mut [usize],
) {
    assert_eq!(data.len(), assignments.len());
    data.iter().zip(assignments.iter_mut()).for_each(|(xs, c)| {
        let (new_center, _) = centers
            .iter()
            .map(|cs| dist(xs.borrow(), cs.borrow()))
            .enumerate()
            .min_by(|x, y| (x.1.partial_cmp(&(y.1)).unwrap()))
            .unwrap();
        *c = new_center;
    });
}
fn update_centers<D: std::borrow::Borrow<[f64]>>(
    data: &[D],
    centers: &mut [Vec<f64>],
    counts: &mut [usize],
    assignments: &[usize],
) {
    centers
        .iter_mut()
        .for_each(|cs| cs.iter_mut().for_each(|x| *x = 0f64));
    counts.iter_mut().for_each(|c| *c = 0);
    for (&asn, xs) in assignments.iter().zip(data.iter()) {
        let xs = xs.borrow();
        centers[asn].iter_mut().zip(xs).for_each(|(c, x)| *c += x);
        counts[asn] += 1;
    }
    centers
        .iter_mut()
        .zip(counts.iter())
        .filter(|&(_, &cou)| 0 < cou)
        .for_each(|(cen, &cou)| cen.iter_mut().for_each(|x| *x /= cou as f64));
}
fn get_dist<D: std::borrow::Borrow<[f64]>>(
    data: &[D],
    centers: &[Vec<f64>],
    assignments: &[usize],
) -> f64 {
    data.iter()
        .zip(assignments)
        .map(|(xs, &asn)| dist(xs.borrow(), &centers[asn]))
        .sum()
}
fn dist(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(xs.len(), ys.len());
    std::iter::zip(xs.iter(), ys.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum()
}

fn suggest_first<R: Rng, D: std::borrow::Borrow<[f64]>>(
    data: &[D],
    k: usize,
    rng: &mut R,
) -> Vec<usize> {
    use rand::seq::SliceRandom;
    assert!(k <= data.len());
    let mut centers: Vec<&[f64]> = vec![data.choose(rng).unwrap().borrow()];
    let choices: Vec<_> = (0..data.len()).collect();
    for _ in 0..k - 1 {
        let dists: Vec<_> = data
            .iter()
            .map(|xs| {
                centers
                    .iter()
                    .map(|cs| dist(xs.borrow(), cs))
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap()
            })
            .collect();
        let idx = *choices.choose_weighted(rng, |&i| dists[i]).unwrap();
        centers.push(data[idx].borrow());
    }
    let mut assignments = vec![0; data.len()];
    update_assignments(data, &centers, &mut assignments);
    assignments
}

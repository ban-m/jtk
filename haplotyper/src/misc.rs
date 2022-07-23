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

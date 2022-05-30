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

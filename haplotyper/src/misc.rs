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

pub fn adjusted_rand_index(label: &[usize], pred: &[usize]) -> f64 {
    assert_eq!(label.len(), pred.len());
    let lab_max = *label.iter().max().unwrap();
    let pred_max = *pred.iter().max().unwrap();
    let mut cont_table = vec![vec![0; pred_max + 1]; lab_max + 1];
    let mut lab_sum = vec![0; lab_max + 1];
    let mut pred_sum = vec![0; pred_max + 1];
    for (&lab, &pred) in label.iter().zip(pred.iter()) {
        cont_table[lab][pred] += 1;
        lab_sum[lab] += 1;
        pred_sum[pred] += 1;
    }
    fn choose(x: &usize) -> usize {
        (x.max(&1) - 1) * x / 2
    }
    let lab_match: usize = lab_sum.iter().map(choose).sum();
    let pred_match: usize = pred_sum.iter().map(choose).sum();
    let num_of_pairs = choose(&label.len());
    let both_match: usize = cont_table.iter().flatten().map(choose).sum();
    assert!(both_match <= (lab_match + pred_match) / 2);
    let match_prod = (lab_match * pred_match) as i64;
    let denom = (num_of_pairs * (lab_match + pred_match) / 2) as i64 - match_prod;
    let numer = (num_of_pairs * both_match) as i64 - match_prod;
    numer as f64 / denom as f64
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

/// Return the k-mer entropy, i.e., \sum_{kmer} -1 *  frac_of_kmer * ln (frac_of_kmer)
pub fn entropy(seq: &[u8], k: usize) -> f64 {
    if seq.len() < k {
        return 0f64;
    } else {
        use std::collections::HashMap;
        let mut counts: HashMap<_, u32> = HashMap::new();
        for kmer in seq.windows(k) {
            let id = crate::repeat_masking::to_idx(kmer);
            *counts.entry(id).or_default() += 1;
        }
        let total = counts.values().sum::<u32>() as f64;
        counts
            .values()
            .map(|&count| {
                assert!(count > 0);
                let f = (count as f64) / total;
                f * f.ln() * -1f64
            })
            .sum()
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

pub fn kiley_op_to_ops(k_ops: &[kiley::Op]) -> definitions::Ops {
    use definitions::Op;
    use definitions::Ops;
    if k_ops.is_empty() {
        return Ops(Vec::new());
    }
    fn is_the_same(op1: kiley::Op, op2: kiley::Op) -> bool {
        match (op1, op2) {
            (kiley::Op::Match, kiley::Op::Mismatch) | (kiley::Op::Mismatch, kiley::Op::Match) => {
                true
            }
            (x, y) if x == y => true,
            _ => false,
        }
    }
    assert!(!k_ops.is_empty());
    let (mut current_op, mut len) = (k_ops[0], 1);
    let mut ops = vec![];
    for &op in k_ops.iter().skip(1) {
        if is_the_same(op, current_op) {
            len += 1;
        } else {
            match current_op {
                kiley::Op::Del => ops.push(Op::Del(len)),
                kiley::Op::Ins => ops.push(Op::Ins(len)),
                kiley::Op::Mismatch | kiley::Op::Match => ops.push(Op::Match(len)),
            }
            current_op = op;
            len = 1;
        }
    }
    match current_op {
        kiley::Op::Del => ops.push(Op::Del(len)),
        kiley::Op::Ins => ops.push(Op::Ins(len)),
        kiley::Op::Mismatch | kiley::Op::Match => ops.push(Op::Match(len)),
    }
    Ops(ops)
}

use definitions::DataSet;
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
        assert!(new_dist < dist + UPDATE_THR, "{},{}", dist, new_dist);
        if dist - new_dist < UPDATE_THR {
            break;
        } else {
            dist = new_dist;
        }
    }
    // for (i, (center, c)) in centers.iter().zip(counts.iter()).enumerate() {
    //     let cs: Vec<_> = center.iter().map(|x| format!("{x:.1}")).collect();
    //     debug!("KMEANS\t{c}\t{i}\t{}", cs.join("\t"));
    // }
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

#[derive(Debug, Clone, PartialEq, Eq)]
struct Fix {
    position: usize,
    base: u8,
    length: usize,
    op_type: kiley::Op,
}

impl Fix {
    // fn hit(&self, (qpos, rpos): (usize, usize)) -> bool {
    //     match self.op_type {
    //         kiley::Op::Ins => qpos == self.position,
    //         kiley::Op::Del => rpos == self.position,
    //         _ => unreachable!(),
    //     }
    // }
    fn new(position: usize, base: u8, length: usize, op_type: kiley::Op) -> Self {
        Self {
            position,
            base,
            length,
            op_type,
        }
    }
}
#[derive(Debug, Clone)]
pub struct Fixes(Vec<Fix>);

/// Fix long homopolymer insertion/deletion.
/// Return the information about the location the truncation/elongation happened.
pub fn fix_long_homopolymers(
    query: &mut Vec<u8>,
    refr: &[u8],
    ops: &mut Vec<kiley::Op>,
    max_allowed: usize,
) -> Fixes {
    use kiley::Op;
    let orig_q_len = query.len();
    let orig_op_len = ops.len();
    let (mut qpos, mut rpos) = (0, 0);
    let mut state: Option<(Op, u8, usize)> = None;
    let mut fixes = vec![];
    for op_pos in 0..orig_op_len {
        state = match (state, ops[op_pos]) {
            (Some((Op::Del, base, len)), op @ (Op::Match | Op::Mismatch)) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - op_len;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Del))
                }
                ops.push(op);
                query.push(query[qpos]);
                rpos += 1;
                qpos += 1;
                None
            }
            (Some((Op::Del, base, len)), Op::Del) if base == refr[rpos] => {
                if len < max_allowed {
                    ops.push(Op::Del);
                } else {
                    ops.push(Op::Match);
                    query.push(base);
                }
                rpos += 1;
                Some((Op::Del, base, len + 1))
            }
            (Some((Op::Del, base, len)), Op::Del) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - op_len;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Del))
                }
                rpos += 1;
                ops.push(Op::Del);
                Some((Op::Del, refr[rpos - 1], 1))
            }
            (Some((Op::Del, base, len)), Op::Ins) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - op_len;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Del))
                }
                let qbase = query[qpos];
                ops.push(Op::Ins);
                query.push(qbase);
                qpos += 1;
                Some((Op::Ins, qbase, 1))
            }
            (Some((Op::Ins, base, len)), op @ (Op::Match | Op::Mismatch)) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - max_allowed;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Ins));
                }
                ops.push(op);
                query.push(query[qpos]);
                qpos += 1;
                rpos += 1;
                None
            }
            (Some((Op::Ins, base, len)), Op::Ins) if base == query[qpos] => {
                if len < max_allowed {
                    ops.push(Op::Ins);
                    query.push(query[qpos]);
                }
                qpos += 1;
                Some((Op::Ins, base, len + 1))
            }
            (Some((Op::Ins, base, len)), Op::Ins) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - max_allowed;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Ins));
                }
                let qbase = query[qpos];
                ops.push(Op::Ins);
                query.push(qbase);
                qpos += 1;
                Some((Op::Ins, qbase, 1))
            }
            (Some((Op::Ins, base, len)), Op::Del) => {
                if max_allowed < len {
                    let op_len = len - max_allowed;
                    let start_pos = query.len() - orig_q_len - max_allowed;
                    fixes.push(Fix::new(start_pos, base, op_len, Op::Ins));
                }
                ops.push(Op::Del);
                rpos += 1;
                Some((Op::Del, refr[rpos - 1], 1))
            }
            (None, op @ (Op::Match | Op::Mismatch)) => {
                ops.push(op);
                query.push(query[qpos]);
                qpos += 1;
                rpos += 1;
                None
            }
            (None, Op::Del) => {
                ops.push(Op::Del);
                rpos += 1;
                Some((Op::Del, refr[rpos - 1], 1))
            }
            (None, Op::Ins) => {
                ops.push(Op::Ins);
                let qbase = query[qpos];
                query.push(qbase);
                qpos += 1;
                Some((Op::Ins, qbase, 1))
            }
            (Some((Op::Mismatch | Op::Match, _, _)), _) => unreachable!(),
        };
    }
    assert_eq!(orig_q_len, qpos);
    assert_eq!(refr.len(), rpos);
    match state {
        Some((Op::Ins, base, len)) if max_allowed < len => {
            let op_len = len - max_allowed;
            let start_pos = query.len() - orig_q_len - op_len;
            fixes.push(Fix::new(start_pos, base, op_len, Op::Ins));
        }
        Some((Op::Del, base, len)) if max_allowed < len => {
            let op_len = len - max_allowed;
            let start_pos = query.len() - orig_q_len - op_len;
            fixes.push(Fix::new(start_pos, base, op_len, Op::Del))
        }
        _ => {}
    }
    // Remove original sequences.
    {
        let mut idx = 0;
        query.retain(|_| {
            idx += 1;
            orig_q_len < idx
        });
        let mut idx = 0;
        ops.retain(|_| {
            idx += 1;
            orig_op_len < idx
        });
    }
    // Check consistency.
    let reflen = ops.iter().filter(|&&op| op != Op::Ins).count();
    let querylen = ops.iter().filter(|&&op| op != Op::Del).count();
    assert_eq!(reflen, refr.len());
    assert_eq!(querylen, query.len());
    Fixes(fixes)
}

/// Recover the original sequence.
pub fn revert_fix_homopolymers(
    query: &mut Vec<u8>,
    ops: &mut Vec<kiley::Op>,
    Fixes(ref modifs): &Fixes,
) {
    let refrlen = ops.iter().filter(|&&op| op != Op::Ins).count();
    let orig_q_len = query.len();
    let orig_op_len = ops.len();
    use kiley::Op;
    let mut qpos = 0;
    let mut skip = 0;
    let mut modifs = modifs.iter();
    let mut next_modif = modifs.next();
    for op_pos in 0..orig_op_len {
        match next_modif {
            Some(fix) if fix.position == qpos => {
                match fix.op_type {
                    Op::Del => {
                        skip += fix.length;
                    }
                    Op::Ins => {
                        ops.extend(std::iter::repeat(Op::Ins).take(fix.length));
                        query.extend(std::iter::repeat(fix.base).take(fix.length));
                    }
                    _ => unreachable!(),
                }
                next_modif = modifs.next();
            }
            _ => {}
        }
        match ops[op_pos] {
            Op::Mismatch | Op::Match if 0 < skip => {
                qpos += 1;
                skip -= 1;
                ops.push(Op::Del);
            }
            Op::Ins if 0 < skip => {
                skip -= 1;
                qpos += 1;
            }
            Op::Mismatch | Op::Match | Op::Ins => {
                query.push(query[qpos]);
                ops.push(ops[op_pos]);
                qpos += 1;
            }
            Op::Del => {
                ops.push(ops[op_pos]);
            }
        }
    }
    // Fix remaining operations.
    loop {
        match next_modif {
            Some(fix) if fix.position == qpos => {
                match fix.op_type {
                    Op::Del => {
                        skip += fix.length;
                    }
                    Op::Ins => {
                        ops.extend(std::iter::repeat(Op::Ins).take(fix.length));
                        query.extend(std::iter::repeat(fix.base).take(fix.length));
                    }
                    _ => unreachable!(),
                }
                next_modif = modifs.next();
            }
            _ => break,
        }
    }
    // Revemo original sequences.
    {
        let mut idx = 0;
        query.retain(|_| {
            idx += 1;
            orig_q_len < idx
        });
        let mut idx = 0;
        ops.retain(|_| {
            idx += 1;
            orig_op_len < idx
        });
    }
    let reflen = ops.iter().filter(|&&op| op != Op::Ins).count();
    let querylen = ops.iter().filter(|&&op| op != Op::Del).count();
    assert_eq!(reflen, refrlen);
    assert_eq!(querylen, query.len());
}

pub fn check_compress_and_recover(ds: &DataSet) {
    use std::collections::HashMap;
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        let unit = units.get(&node.unit).unwrap();
        let (_, ar, _) = node.recover(unit);
        let orig_iden = ar.iter().filter(|&&x| x == b'|').count() as f64 / ar.len() as f64;
        let orig_seq = node.seq().to_vec();
        let orig_ops = ops_to_kiley(&node.cigar);
        let mut ops = orig_ops.clone();
        let mut seq = orig_seq.clone();
        let fixes = fix_long_homopolymers(&mut seq, unit.seq(), &mut ops, 2);
        revert_fix_homopolymers(&mut seq, &mut ops, &fixes);
        let seq = std::str::from_utf8(&seq).unwrap();
        let orig_seq = std::str::from_utf8(&orig_seq).unwrap();
        if seq != orig_seq {
            let ops = kiley::bialignment::edit_dist_ops(orig_seq.as_bytes(), seq.as_bytes()).1;
            let (xr, ar, yr) = kiley::recover(orig_seq.as_bytes(), seq.as_bytes(), &ops);
            for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
                eprintln!("{}", std::str::from_utf8(xr).unwrap());
                eprintln!("{}", std::str::from_utf8(ar).unwrap());
                eprintln!("{}\n", std::str::from_utf8(yr).unwrap());
            }
        }
        assert_eq!(seq, orig_seq);
        let (_, ar, _) = node.recover(unit);
        let iden = ar.iter().filter(|&&x| x == b'|').count() as f64 / ar.len() as f64;
        assert!(orig_iden - iden < 0.05);
    }
}

// The maximum value of sum of a range in xs,
// If the sequence is empty, return i64::MIN
pub fn max_region<T: std::iter::Iterator<Item = i64>>(xs: T) -> i64 {
    // The max value of the sum of the range ending at i.
    let mut right = i64::MIN;
    // The max value of the sum of the range ending before i.
    let mut left = i64::MIN;
    for x in xs {
        left = right.max(left);
        right = match right.is_negative() {
            true => x,
            false => right + x,
        };
    }
    right.max(left)
}

/// Return the max_{i<j} xs[i..j], where xs is the converted value of the alignment operation.
/// Each deletion and insertions would be converted into `indel_weight * operation length`, whereas
/// match would be into `mat_weight * operation length * -1`.
/// Thus, if we set mat_weight sufficiently large and indel_weight = 1, this function would return
/// the maximum length of the consective insertion/deletion.
pub fn max_indel(ops: &definitions::Ops, mat_weight: u64, indel_weight: u64) -> u64 {
    let (mat_weight, indel_weight) = (mat_weight as i64, indel_weight as i64);
    let xs = ops.0.iter().map(|&op| match op {
        definitions::Op::Match(l) => -(l as i64) * mat_weight,
        definitions::Op::Del(l) => l as i64 * indel_weight,
        definitions::Op::Ins(l) => l as i64 * indel_weight,
    });
    max_region(xs).max(0) as u64
}

pub fn max_indel_cigar(cigar: &[bio_utils::sam::Op], mat_weight: u64, indel_weight: u64) -> u64 {
    let (mat_weight, indel_weight) = (mat_weight as i64, indel_weight as i64);
    let xs = cigar.iter().map(|&op| match op {
        bio_utils::sam::Op::Align(l)
        | bio_utils::sam::Op::Match(l)
        | bio_utils::sam::Op::Mismatch(l) => -(l as i64) * mat_weight,
        bio_utils::sam::Op::Insertion(l)
        | bio_utils::sam::Op::Deletion(l)
        | bio_utils::sam::Op::SoftClip(l)
        | bio_utils::sam::Op::HardClip(l) => (l as i64) * indel_weight,
        _ => 0,
    });
    max_region(xs).max(0) as u64
}

pub fn max_indel_node(node: &definitions::Node, mat_weight: u64, indel_weight: u64) -> u64 {
    max_indel(&node.cigar, mat_weight, indel_weight)
}

pub fn update_coverage(ds: &mut DataSet) {
    if !ds.coverage.is_protected() {
        use std::collections::HashMap;
        let mut counts: HashMap<_, u32> = HashMap::new();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *counts.entry(node.unit).or_default() += 1;
        }
        let mut counts: Vec<_> = counts.values().copied().collect();
        counts.sort_unstable();
        let cov = counts[counts.len() / 2] as f64 / 2f64;
        debug!("MULTP\tCOVERAGE\t{}\tHAPLOID", cov);
        ds.coverage.set(cov);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kiley::Op;
    #[test]
    fn homop_compress() {
        let mut query = b"AAAAAACCTTTTTGTT".to_vec();
        let refr = b"CCTTCCCCCCCCCGTT";
        let mut ops = vec![
            vec![Op::Ins; 6],
            vec![Op::Match; 4],
            vec![Op::Ins; 3],
            vec![Op::Del; 9],
            vec![Op::Match; 3],
        ]
        .concat();
        let original_q = String::from_utf8(query.clone()).unwrap();
        let original_op = ops.clone();
        let fix_ops = fix_long_homopolymers(&mut query, refr, &mut ops, 2);
        let query_str = std::str::from_utf8(&query).unwrap();
        assert_eq!(query_str, "AACCTTTTCCCCCCCGTT");
        let ans_ops = vec![
            vec![Op::Ins; 2],
            vec![Op::Match; 4],
            vec![Op::Ins; 2],
            vec![Op::Del; 2],
            vec![Op::Match; 7],
            vec![Op::Match; 3],
        ]
        .concat();
        assert_eq!(ans_ops, ops);
        let ans_fix_ops = vec![
            Fix::new(0, b'A', 4, Op::Ins),
            Fix::new(6, b'T', 1, Op::Ins),
            Fix::new(8, b'C', 7, Op::Del),
        ];
        assert_eq!(ans_fix_ops, fix_ops.0);
        revert_fix_homopolymers(&mut query, &mut ops, &fix_ops);
        assert_eq!(ops, original_op);
        assert_eq!(original_q, std::str::from_utf8(&query).unwrap());
    }

    #[test]
    fn homop_compress2() {
        let mut query = b"GTAATAAAAAA".to_vec();
        let refr = b"GCACAGAA";
        let mut ops = vec![
            vec![Op::Match, Op::Mismatch, Op::Match],
            vec![Op::Del],
            vec![Op::Ins; 4],
            vec![Op::Match, Op::Mismatch, Op::Match, Op::Match],
        ]
        .concat();
        let original_q = String::from_utf8(query.clone()).unwrap();
        let original_op = ops.clone();
        let fix_ops = fix_long_homopolymers(&mut query, refr, &mut ops, 2);
        eprintln!("{:?}", fix_ops.0);
        revert_fix_homopolymers(&mut query, &mut ops, &fix_ops);
        assert_eq!(ops, original_op);
        assert_eq!(original_q, std::str::from_utf8(&query).unwrap());
    }
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    #[test]
    fn homop_compress_check_random() {
        let profile = kiley::gen_seq::Profile::new(0.05, 0.05, 0.05);
        for seed in 0..100 {
            // println!("{seed}");
            let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
            let seq = kiley::gen_seq::generate_seq(&mut rng, 100);
            let query = kiley::gen_seq::introduce_randomness(&seq, &mut rng, &profile);
            // println!("{}", std::str::from_utf8(&seq).unwrap());
            // println!("{}", std::str::from_utf8(&query).unwrap());
            let ops = kiley::bialignment::edit_dist_ops(&seq, &query).1;
            let mut modif_q = query.clone();
            let mut modif_ops = ops.clone();
            let fix_ops = fix_long_homopolymers(&mut modif_q, &seq, &mut modif_ops, 2);
            revert_fix_homopolymers(&mut modif_q, &mut modif_ops, &fix_ops);
            assert_eq!(ops, modif_ops);
            let query = std::str::from_utf8(&query).unwrap();
            let modif_q = std::str::from_utf8(&modif_q).unwrap();
            assert_eq!(query, modif_q);
        }
    }
    #[test]
    fn max_range_operation_test() {
        use definitions::Op;
        let ops = [
            Op::Match(10),
            Op::Del(5),
            Op::Ins(2),
            Op::Match(1),
            Op::Del(5),
            Op::Match(10),
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Ins(l) => -(*l as i64),
            Op::Del(l) => *l as i64 * 2,
        });
        let max_del = max_region(iter);
        assert_eq!(max_del, 10 + 10 - 3);
        let ops = [
            Op::Ins(10),
            Op::Del(5),
            Op::Ins(2),
            Op::Match(1),
            Op::Del(5),
            Op::Match(10),
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Del(l) => -(*l as i64),
            Op::Ins(l) => *l as i64 * 2,
        });
        let max_in = max_region(iter);
        assert_eq!(max_in, 20);
        let ops = [
            Op::Ins(10),
            Op::Del(5),
            Op::Ins(2),    // 19
            Op::Match(1),  // 18
            Op::Del(5),    // 13
            Op::Match(10), // 3
            Op::Ins(100),  // 203
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Del(l) => -(*l as i64),
            Op::Ins(l) => *l as i64 * 2,
        });
        let max_in = max_region(iter);
        assert_eq!(max_in, 203);
    }
    #[test]
    fn rand_index_test() {
        let pred = [0, 0, 0, 1, 1, 1];
        let answ = [0, 0, 1, 1, 2, 2];
        assert!((0.6666 - rand_index(&pred, &answ)).abs() < 0.0001);
    }
}

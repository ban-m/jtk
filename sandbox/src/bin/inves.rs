use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
const MISMATCH: f64 = 0.15;
const MATCH: f64 = 1f64 - MISMATCH;
const MISM_SCORE: f64 = -1.897;
const GAP_SCORE: f64 = -4.605;
const MISM_UNIT: f64 = -10000000f64;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    eprintln!("{:?}", std::time::Instant::now() - start);
    let alignment_score = get_match_score(&ds);
    for ((unit, cl), score) in alignment_score.iter() {
        eprintln!("SCORE\t{}\t{}\t{:.3}", unit, cl, score);
    }
    eprintln!("ALN\tLEN1\tLEN2\tScore\tNaive\tMM");
    for (i, read1) in ds.encoded_reads.iter().enumerate() {
        let read1: Vec<_> = read1.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        for read2 in ds.encoded_reads.iter().skip(i + 1) {
            let mut read2: Vec<_> = read2.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
            let forward = aln_score_ops(&read1, &read2, &alignment_score);
            read2.reverse();
            let reverse = aln_score_ops(&read1, &read2, &alignment_score);
            let ((score, ops), _is_forward) = if reverse.0 < forward.0 {
                (forward, true)
            } else {
                (reverse, false)
            };
            let reverse = edit_dist(&read1, &read2);
            read2.reverse();
            let forward = edit_dist(&read1, &read2);
            let edit = reverse.min(forward);
            if 0f64 < score {
                let mm = ops
                    .iter()
                    .filter(|x| match x {
                        Cigar::Match(_, _) => true,
                        _ => false,
                    })
                    .count();
                eprintln!(
                    "ALN\t{}\t{}\t{}\t{}\t{}",
                    read1.len(),
                    read2.len(),
                    score,
                    edit,
                    mm
                );
            }
            // println!("{:?}\n{:?}\n{:?}\t{}\n", read1, read2, ops, is_forward);
        }
    }
    Ok(())
}

fn edit_dist(xs: &[(u64, u64)], ys: &[(u64, u64)]) -> i32 {
    let mut dp = vec![vec![0; ys.len() + 1]; xs.len() + 1];
    for (i, x) in xs.iter().enumerate().map(|(i, p)| (i + 1, p)) {
        for (j, y) in ys.iter().enumerate().map(|(i, p)| (i + 1, p)) {
            let mat = if x == y { 1 } else { -1 };
            dp[i][j] = (dp[i - 1][j - 1] + mat)
                .max(dp[i - 1][j] - 1)
                .max(dp[i][j - 1] - 1);
        }
    }
    (0..xs.len() + 1)
        .map(|i| (i, ys.len()))
        .chain((0..ys.len()).map(|j| (xs.len(), j)))
        .map(|(i, j)| dp[i][j])
        .max()
        .unwrap()
}

#[derive(Clone, PartialEq, Eq)]
enum Cigar {
    Match(u64, u64),
    Ins(u64, u64),
    Del,
}

impl std::fmt::Debug for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Cigar::Match(_, _) => write!(f, "M"),
            Cigar::Ins(_, _) => write!(f, "I"),
            Cigar::Del => write!(f, "D"),
        }
    }
}

fn aln_score_ops(
    xs: &[(u64, u64)],
    ys: &[(u64, u64)],
    score: &HashMap<(u64, u64), f64>,
) -> (f64, Vec<Cigar>) {
    let mut dp = vec![vec![0f64; ys.len() + 1]; xs.len() + 1];
    for (i, (u1, c1)) in xs.iter().enumerate().map(|(i, &p)| (i + 1, p)) {
        for (j, (u2, c2)) in ys.iter().enumerate().map(|(j, &p)| (j + 1, p)) {
            let mat_score = match (u1 == u2, c1 == c2) {
                (true, true) => score[&(u1, c1)],
                (true, false) => MISM_SCORE,
                (false, _) => MISM_UNIT,
            };
            dp[i][j] = (dp[i - 1][j - 1] + mat_score)
                .max(dp[i - 1][j] + GAP_SCORE)
                .max(dp[i][j - 1] + GAP_SCORE);
        }
    }
    let (opt_score, mut i, mut j) = (0..xs.len() + 1)
        .map(|i| (i, ys.len()))
        .chain((0..ys.len()).map(|j| (xs.len(), j)))
        .map(|(i, j)| (dp[i][j], i, j))
        .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
        .unwrap();
    let mut ops = vec![];
    for jpos in (j..ys.len()).rev() {
        let (u2, c2) = ys[jpos];
        ops.push(Cigar::Ins(u2, c2));
    }
    for _ in (i..xs.len()).rev() {
        ops.push(Cigar::Del);
    }
    while 0 < i && 0 < j {
        let (u1, c1) = xs[i - 1];
        let (u2, c2) = ys[j - 1];
        let mat_score = match (u1 == u2, c1 == c2) {
            (true, true) => score[&(u1, c1)],
            (true, false) => MISM_SCORE,
            (false, _) => MISM_UNIT,
        };
        if (dp[i][j] - (dp[i - 1][j - 1] + mat_score)).abs() < 0.00001 {
            ops.push(Cigar::Match(u2, c2));
            i -= 1;
            j -= 1;
        } else if (dp[i][j] - (dp[i - 1][j] + GAP_SCORE)).abs() < 0.00001 {
            ops.push(Cigar::Del);
            i -= 1;
        } else if (dp[i][j] - (dp[i][j - 1] + GAP_SCORE)).abs() < 0.00001 {
            ops.push(Cigar::Ins(u2, c2));
            j -= 1;
        } else {
            unreachable!()
        }
    }
    for jpos in (0..j).rev() {
        let (u2, c2) = ys[jpos];
        ops.push(Cigar::Ins(u2, c2));
    }
    ops.extend(std::iter::repeat(Cigar::Del).take(i));
    ops.reverse();
    (opt_score, ops)
}

fn get_match_score(ds: &DataSet) -> HashMap<(u64, u64), f64> {
    let mut counts: HashMap<u64, HashMap<u64, u32>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *counts
                .entry(node.unit)
                .or_default()
                .entry(node.cluster)
                .or_default() += 1;
        }
    }
    counts
        .iter()
        .flat_map(|(&unit, val)| {
            let total: u32 = val.values().sum();
            val.iter()
                .map(|(&cluster, &count)| {
                    let frac = count as f64 / total as f64;
                    let score = (MATCH + MISMATCH * frac).ln() - frac.ln();
                    ((unit, cluster), score)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

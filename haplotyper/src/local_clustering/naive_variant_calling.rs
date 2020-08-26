use super::ChunkedUnit;
use super::ClusteringConfig;
use poa_hmm::POA;
use rand::Rng;
// use rayon::prelude::*;

fn select_variants(
    mut variants: Vec<Vec<Vec<f64>>>,
    chain_len: usize,
    variant_number: usize,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>) {
    let mut position = vec![false; chain_len];
    let thr = {
        let mut var: Vec<_> = variants
            .iter()
            .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
            .copied()
            .collect();
        var.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let pos = var.len() - variant_number;
        var[pos].max(0.05)
    };
    for bss in variants.iter_mut() {
        for bs in bss.iter_mut() {
            for (idx, b) in bs.iter_mut().enumerate() {
                if *b < thr {
                    *b = 0.;
                } else {
                    position[idx] = true;
                }
            }
        }
    }
    (variants, position)
}

pub fn get_variants_alignment<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>, f64) {
    let betas = variant_calling_by_alignment(data, chain_len, rng, c);
    // for bss in betas.iter() {
    //     for bs in bss.iter() {
    //         let line: Vec<_> = bs
    //             .iter()
    //             .enumerate()
    //             .map(|(i, b)| format!("{}:{:.3}", i, b))
    //             .collect();
    //         debug!("{:?}", line.join(","));
    //     }
    // }
    let (betas, position_in_use) = select_variants(betas, chain_len, c.variant_num);
    // let pos: Vec<_> = position_in_use
    //     .iter()
    //     .enumerate()
    //     .filter(|&(_, &b)| b)
    //     .map(|x| x.0)
    //     .collect();
    // debug!("{:?}", pos);
    (betas, position_in_use, 0.)
}

pub fn variant_calling_by_alignment<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> Vec<Vec<Vec<f64>>> {
    (0..c.cluster_num)
        .map(|i| {
            (0..i)
                .map(|j| {
                    let mut pileup = vec![vec![]; chain_len];
                    for read in data.iter().filter(|d| d.cluster == i || d.cluster == j) {
                        for chunk in read.chunks.iter() {
                            pileup[chunk.pos].push(chunk.seq.as_slice());
                        }
                    }
                    use rand::seq::SliceRandom;
                    pileup.iter_mut().for_each(|seqs| seqs.shuffle(rng));
                    let mut betas: Vec<_> =
                        pileup.iter().map(variant_calling_from_pileup).collect();
                    let sum = betas.iter().map(|x| x * x).sum::<f64>().sqrt();
                    if sum < 0.0001 {
                        vec![(pileup.len() as f64).recip(); pileup.len()]
                    } else {
                        betas.iter_mut().for_each(|x| *x /= sum);
                        betas
                    }
                })
                .collect()
        })
        .collect()
}

fn variant_calling_from_pileup(seqs: &Vec<&[u8]>) -> f64 {
    let consensus = POA::from_slice_default(&seqs).consensus();
    let pileup = seqs
        .iter()
        .filter_map(|seq| alignment(seq, &consensus, (1, -1, -1)))
        .fold(Pileup::new(&consensus), |x, aln| x.add(aln.1));
    pileup.column.iter().map(|p| p.variant_call()).sum::<f64>()
}

fn score(x: u8, y: u8, mat: i32, mism: i32) -> i32 {
    if x == y {
        mat
    } else {
        mism
    }
}

fn alignment(
    qry: &[u8],
    rfr: &[u8],
    (mat, mism, gap): (i32, i32, i32),
) -> Option<(i32, Vec<Cigar>)> {
    let mut dp = vec![vec![0; rfr.len() + 1]; qry.len() + 1];
    for (i, &q) in qry.iter().enumerate() {
        for (j, &r) in rfr.iter().enumerate() {
            let mat = dp[i][j] + score(q, r, mat, mism);
            let ins = dp[i][j + 1] + gap;
            let del = dp[i + 1][j] + gap;
            dp[i + 1][j + 1] = mat.max(ins).max(del);
        }
    }
    // Determine the starting point.
    let (row_pos, row_max) = dp.last()?.iter().enumerate().max_by_key(|x| x.1)?;
    let (column_pos, column_max) = dp
        .iter()
        .filter_map(|x| x.last())
        .enumerate()
        .max_by_key(|x| x.1)?;
    let score = *column_max.max(row_max);
    if score <= 0 {
        return None;
    }
    let (mut q_pos, mut r_pos) = if row_max < column_max {
        (column_pos, rfr.len())
    } else {
        (qry.len(), row_pos)
    };
    assert_eq!(dp[q_pos][r_pos], *column_max.max(row_max));
    let mut cigar = vec![];
    for q in (q_pos..qry.len()).rev() {
        cigar.push(Cigar::Ins(qry[q]));
    }
    for _ in (r_pos..rfr.len()).rev() {
        cigar.push(Cigar::Del);
    }
    // Traceback.
    while q_pos > 0 && r_pos > 0 {
        let current = dp[q_pos][r_pos];
        let op = if current == dp[q_pos - 1][r_pos] + gap {
            let base = qry[q_pos - 1];
            q_pos -= 1;
            Cigar::Ins(base)
        } else if current == dp[q_pos][r_pos - 1] + gap {
            r_pos -= 1;
            Cigar::Del
        } else {
            let base = qry[q_pos - 1];
            r_pos -= 1;
            q_pos -= 1;
            Cigar::Match(base)
        };
        cigar.push(op);
    }
    while q_pos > 0 {
        let base = qry[q_pos - 1];
        cigar.push(Cigar::Ins(base));
        q_pos -= 1;
    }
    while r_pos > 0 {
        cigar.push(Cigar::Del);
        r_pos -= 1;
    }
    cigar.reverse();
    Some((score, cigar))
}

#[derive(Clone, PartialEq, Eq)]
enum Cigar {
    Match(u8),
    Ins(u8),
    Del,
}

impl std::fmt::Debug for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Cigar::Match(_) => write!(f, "M"),
            Cigar::Ins(_) => write!(f, "I"),
            Cigar::Del => write!(f, "D"),
        }
    }
}

#[derive(Clone, Debug)]
struct Pileup {
    column: Vec<Column>,
}

impl Pileup {
    fn new(ns: &[u8]) -> Self {
        let column: Vec<_> = ns.iter().map(|&n| Column::new(n)).collect();
        Self { column }
    }
    fn add(mut self, y: Vec<Cigar>) -> Self {
        let q_len = y
            .iter()
            .map(|op| match op {
                Cigar::Match(_) => 1,
                Cigar::Ins(_) => 1,
                _ => 0,
            })
            .sum::<usize>();
        let (mut r_pos, mut q_pos) = (0, 0);
        for op in y {
            match op {
                Cigar::Match(base) => {
                    self.column[r_pos].m[b2i(base)] += 1;
                    self.column[r_pos].c += 1;
                    r_pos += 1;
                    q_pos += 1;
                }
                Cigar::Ins(base) => {
                    if 0 < r_pos && r_pos < self.column.len() {
                        self.column[r_pos].i.push(base);
                    }
                    q_pos += 1;
                }
                Cigar::Del => {
                    if 0 < q_pos && q_pos < q_len {
                        self.column[r_pos].d += 1;
                    }
                    r_pos += 1;
                }
            }
        }
        self
    }
}

#[derive(Clone, Debug)]
struct Column {
    m: [u8; 4],
    i: Vec<u8>,
    d: usize,
    c: usize,
}

fn b2i(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

impl Column {
    fn new(b: u8) -> Self {
        let mut bases = [0; 4];
        bases[b2i(b)] += 1;
        Self {
            m: bases,
            i: vec![],
            d: 0,
            c: 1,
        }
    }
    fn variant_call(&self) -> f64 {
        let max = *self.m.iter().max().unwrap();
        let maf = *self.m.iter().filter(|&&x| x != max).max().unwrap_or(&0);
        let coverage = (self.c + self.i.len()) as f64;
        if maf as f64 / coverage > 0.2 {
            maf as f64 / coverage
        } else {
            0.
        }
        // if maf as f64 / coverage > 0.2 {
        //     maf as f64 / coverage
        // } else if self.i.len() as f64 / coverage > 0.5 {
        //     self.i.len() as f64 / coverage
        // } else if self.d as f64 / coverage > 0.5 {
        //     self.d as f64 / coverage
        // } else {
        //     0.
        // }
    }
}

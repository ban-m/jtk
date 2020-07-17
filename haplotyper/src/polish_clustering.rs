use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct PolishClusteringConfig {
    pub threads: usize,
    pub mat_score: i32,
    pub mismat_score: i32,
    pub gap_score: i32,
}

impl PolishClusteringConfig {
    pub fn new(threads: usize, mat_score: i32, mismat_score: i32, gap_score: i32) -> Self {
        Self {
            threads,
            mat_score,
            mismat_score,
            gap_score,
        }
    }
}
pub trait PolishClustering {
    fn polish_clustering(self, c: &PolishClusteringConfig) -> Self;
}

impl PolishClustering for DataSet {
    fn polish_clustering(mut self, c: &PolishClusteringConfig) -> Self {
        rayon::ThreadPoolBuilder::new()
            .num_threads(c.threads)
            .build_global()
            .unwrap();
        let id_to_cluster: HashMap<_, _> = self
            .assignments
            .iter()
            .map(|asn| (asn.id, asn.cluster))
            .collect();
        let mut buckets: HashMap<_, Vec<&mut EncodedRead>> = HashMap::new();
        for read in self.encoded_reads.iter_mut() {
            if let Some(cluster) = id_to_cluster.get(&read.id) {
                buckets.entry(cluster).or_default().push(read);
            }
        }
        buckets
            .into_iter()
            .for_each(|(_, bucket)| correct(bucket, c));
        self
    }
}

fn score((q_u, q_c): (u64, u64), (r_u, r_c): (u64, u64), mat: i32, mism: i32) -> i32 {
    if q_u == r_u && q_c == r_c {
        mat
    } else if q_u == r_u {
        mism
    } else {
        -100
    }
}

// Align the query to the reference and
// return the edit operations. Note that
// A Cigar::Match(x,y) mean the query sequence at that point is (x,y)
// And Cigar::Ins is a insertion to the reference.
// Also, the alignment is "semi-global" one. See the initialization step.
// TODO: faster!
fn alignment(
    qry: &[(u64, u64)],
    rfr: &[(u64, u64)],
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
    if score <= mat {
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
        let (unit, cluster) = qry[q];
        cigar.push(Cigar::Ins(unit, cluster));
    }
    for _ in (r_pos..rfr.len()).rev() {
        cigar.push(Cigar::Del);
    }
    // Traceback.
    while q_pos > 0 && r_pos > 0 {
        let current = dp[q_pos][r_pos];
        let op = if current == dp[q_pos - 1][r_pos] + gap {
            let (unit, cluster) = qry[q_pos - 1];
            q_pos -= 1;
            Cigar::Ins(unit, cluster)
        } else if current == dp[q_pos][r_pos - 1] + gap {
            r_pos -= 1;
            Cigar::Del
        } else {
            let (unit, cluster) = qry[q_pos - 1];
            r_pos -= 1;
            q_pos -= 1;
            Cigar::Match(unit, cluster)
        };
        cigar.push(op);
    }
    while q_pos > 0 {
        let (unit, cluster) = qry[q_pos - 1];
        cigar.push(Cigar::Ins(unit, cluster));
        q_pos -= 1;
    }
    while r_pos > 0 {
        cigar.push(Cigar::Del);
        r_pos -= 1;
    }
    cigar.reverse();
    Some((score, cigar))
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Cigar {
    Match(u64, u64),
    Ins(u64, u64),
    Del,
}

#[derive(Clone, Debug)]
struct Pileup {
    column: Vec<Column>,
}

impl Pileup {
    fn new(ns: &[(u64, u64)]) -> Self {
        let column: Vec<_> = ns.iter().map(|&n| Column::new(n)).collect();
        Self { column }
    }
    fn add(mut self, y: Vec<Cigar>) -> Self {
        let mut r_pos = 0;
        for op in y {
            match op {
                Cigar::Match(unit, cl) => {
                    self.column[r_pos].m.push((unit, cl));
                    r_pos += 1;
                }
                Cigar::Ins(_, _) => {}
                Cigar::Del => {
                    r_pos += 1;
                }
            }
        }
        self
    }
}

#[derive(Clone, Debug)]
struct Column {
    m: Vec<(u64, u64)>,
}

impl Column {
    fn new(n: (u64, u64)) -> Self {
        Self { m: vec![n] }
    }
    fn generate(&self) -> (u64, u64) {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for &x in self.m.iter() {
            *counts.entry(x).or_default() += 1;
        }
        counts.into_iter().max_by_key(|x| x.1).unwrap().0
    }
}

fn correct(mut reads: Vec<&mut EncodedRead>, config: &PolishClusteringConfig) {
    let reads_summary: Vec<Vec<(u64, u64)>> = reads
        .iter()
        .map(|read| read.nodes.iter().map(|n| (n.unit, n.cluster)).collect())
        .collect();
    let rev_for_reads: Vec<_> = {
        let rev = reads_summary
            .iter()
            .map(|read| read.iter().copied().rev().collect::<Vec<_>>());
        reads_summary.iter().cloned().zip(rev).collect()
    };
    let corrected_reads: Vec<_> = reads_summary
        .par_iter()
        .map(|read| correct_read(read, &rev_for_reads, config))
        .collect();
    assert_eq!(reads.len(), corrected_reads.len());
    for (read, corrected) in reads.iter_mut().zip(corrected_reads) {
        assert_eq!(read.nodes.len(), corrected.len());
        for (node, (unit, cluster)) in read.nodes.iter_mut().zip(corrected) {
            node.unit = unit;
            node.cluster = cluster;
        }
    }
}

fn correct_read(
    read: &Vec<(u64, u64)>,
    reads: &[(Vec<(u64, u64)>, Vec<(u64, u64)>)],
    c: &PolishClusteringConfig,
) -> Vec<(u64, u64)> {
    let param = (c.mat_score, c.mismat_score, c.gap_score);
    let pileup = reads
        .iter()
        .filter_map(|(forward, rev)| match alignment(forward, read, param) {
            Some(res) => Some(res),
            None => alignment(rev, read, param),
        })
        .fold(Pileup::new(read), |x, (_, y)| x.add(y));
    pileup
        .column
        .into_iter()
        .map(|column| column.generate())
        .collect()
}

#[cfg(test)]
mod tests {
    #[derive(Clone, Copy, Debug)]
    struct TestConfig {
        cl: usize,
        num: usize,
        fail: f64,
        max_len: usize,
        min_len: usize,
        unit_len: usize,
    }
    use super::*;
    use rand::Rng;
    use rand_xoshiro::Xoshiro256Plus;
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> Vec<Vec<(u64, u64)>> {
        let TestConfig {
            cl,
            num,
            fail,
            max_len,
            min_len,
            unit_len,
        } = conf;
        (0..num)
            .map(|i| {
                let cluster = (i % cl) as u64;
                let len = r.gen::<usize>() % (max_len - min_len) + min_len;
                let start = r.gen::<usize>() % (unit_len - len);
                let units: Vec<_> = if r.gen_bool(0.5) {
                    (start..=start + len).collect()
                } else {
                    let start = start + len;
                    (start - len..=start).rev().collect()
                };
                units
                    .iter()
                    .map(|&unit| {
                        if r.gen_bool(fail) {
                            let cluster = r.gen::<u64>() % cl as u64;
                            (unit as u64, cluster)
                        } else {
                            (unit as u64, cluster)
                        }
                    })
                    .collect()
            })
            .collect()
    }
    #[test]
    fn error_correction_multi() {
        let mut rng: Xoshiro256Plus = rand::SeedableRng::seed_from_u64(100);
        let conf = TestConfig {
            cl: 2,
            num: 200,
            fail: 0.01,
            max_len: 10,
            min_len: 8,
            unit_len: 40,
        };
        let reads = gen_dataset(&mut rng, conf);
        let rev_for_reads: Vec<_> = {
            let rev = reads
                .iter()
                .map(|read| read.iter().copied().rev().collect::<Vec<_>>());
            reads.iter().cloned().zip(rev).collect()
        };
        let config = PolishClusteringConfig::new(1, 1, -1, -2);
        for read in reads.iter() {
            eprintln!("Correcting:{:?}", read);
            let res = correct_read(&read, &rev_for_reads, &config);
            eprintln!("Corrected:{:?}", res);
            assert_eq!(read.len(), res.len());
            let cl = res[0].1;
            let cluster = res.iter().all(|n| n.1 == cl);
            assert!(cluster, "{:?}", res);
            let suc = res
                .windows(2)
                .all(|w| w[0].0 + 1 == w[1].0 || w[1].0 + 1 == w[0].0);
            assert!(suc, "{:?}", res);
        }
    }
}

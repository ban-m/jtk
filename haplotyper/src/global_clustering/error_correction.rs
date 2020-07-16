use definitions::DataSet;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Clone, Debug)]
pub struct CorrectedRead {
    pub id: u64,
    pub nodes: Vec<Unit>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unit {
    pub unit: u64,
    pub cluster: u64,
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
        let q_len = y
            .iter()
            .map(|op| match op {
                Cigar::Match(_, _) => 1,
                Cigar::Ins(_, _) => 1,
                _ => 0,
            })
            .sum::<usize>();
        let (mut r_pos, mut q_pos) = (0, 0);
        for op in y {
            match op {
                Cigar::Match(unit, cl) => {
                    self.column[r_pos].m.push((unit, cl));
                    self.column[r_pos].c += 1;
                    r_pos += 1;
                    q_pos += 1;
                }
                Cigar::Ins(unit, cl) => {
                    if 0 < r_pos && r_pos < self.column.len() {
                        self.column[r_pos].i.push((unit, cl));
                    }
                    q_pos += 1;
                }
                Cigar::Del => {
                    if 0 < q_pos && q_pos < q_len {
                        self.column[r_pos].c += 1;
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
    m: Vec<(u64, u64)>,
    i: Vec<(u64, u64)>,
    d: usize,
    c: usize,
}

impl Column {
    fn new(n: (u64, u64)) -> Self {
        Self {
            m: vec![n],
            i: vec![],
            d: 0,
            c: 1,
        }
    }
    fn generate(&self, node: &mut Vec<Unit>) {
        if self.i.len() > self.c / 2 {
            let mut counts: HashMap<_, u32> = HashMap::new();
            for &x in self.i.iter() {
                *counts.entry(x).or_default() += 1;
            }
            if let Some(((unit, cluster), _)) = counts.into_iter().max_by_key(|x| x.1) {
                node.push(Unit { unit, cluster });
            }
        }
        if self.m.len() > self.c / 2 {
            let mut counts: HashMap<_, u32> = HashMap::new();
            for &x in self.m.iter() {
                *counts.entry(x).or_default() += 1;
            }
            if let Some(((unit, cluster), _)) = counts.into_iter().max_by_key(|x| x.1) {
                node.push(Unit { unit, cluster });
            }
        }
    }
}
use super::GlobalClusteringConfig;
pub fn local_correction(ds: &DataSet, c: &GlobalClusteringConfig) -> Vec<CorrectedRead> {
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|read| {
            let nodes: Vec<_> = read
                .nodes
                .iter()
                .map(|node| (node.unit, node.cluster))
                .collect();
            (read.id, nodes)
        })
        .collect();
    debug!("Correcting {} reads...", reads.len());
    local_correction_inner(reads, (c.mat_score, c.mismat_score, c.gap_score), 2)
}
fn local_correction_inner(
    reads: Vec<(u64, Vec<(u64, u64)>)>,
    param: (i32, i32, i32),
    thr: i32,
) -> Vec<CorrectedRead> {
    let rev_for: Vec<_> = {
        let mut temp = reads.clone();
        let rev = reads
            .iter()
            .map(|&(id, ref read)| (id, read.iter().rev().copied().collect()));
        temp.extend(rev);
        temp
    };
    reads
        .par_iter()
        .map(|read| correct(read, &rev_for, param, thr))
        .collect()
}

fn correct(
    &(id, ref nodes): &(u64, Vec<(u64, u64)>),
    reads: &[(u64, Vec<(u64, u64)>)],
    param: (i32, i32, i32),
    thr: i32,
) -> CorrectedRead {
    let pileup = reads
        .iter()
        .filter_map(|&(_, ref query)| alignment(query, nodes, param))
        .filter(|&(score, _)| score > param.0 * thr)
        .fold(Pileup::new(nodes), |x, (_, y)| x.add(y));
    let mut nodes = vec![];
    for column in pileup.column {
        column.generate(&mut nodes);
    }
    CorrectedRead { id, nodes }
}

#[cfg(test)]
mod tests {
    #[derive(Clone, Copy, Debug)]
    struct TestConfig {
        cl: usize,
        num: usize,
        fail: f64,
        skip: f64,
        max_len: usize,
        min_len: usize,
        unit_len: usize,
    }
    use super::*;
    use rand::Rng;
    use rand_xoshiro::Xoshiro256Plus;
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> Vec<(u64, Vec<(u64, u64)>)> {
        let TestConfig {
            cl,
            num,
            fail,
            skip,
            max_len,
            min_len,
            unit_len,
        } = conf;
        let mut reads = vec![];
        for i in 0..num {
            let cluster = (i % cl) as u64;
            let len = r.gen::<usize>() % (max_len - min_len) + min_len;
            let start = r.gen::<usize>() % (unit_len - len);
            let units: Vec<_> = if r.gen_bool(0.5) {
                (start..=start + len).collect()
            } else {
                let start = start + len;
                (start - len..=start).rev().collect()
            };
            let mut read = vec![];
            for unit in units {
                if r.gen_bool(skip) {
                    continue;
                } else if r.gen_bool(fail) {
                    let cluster = r.gen::<u64>() % cl as u64;
                    read.push((unit as u64, cluster));
                } else {
                    read.push((unit as u64, cluster));
                }
            }
            reads.push((i as u64, read));
        }
        reads
    }
    #[test]
    fn error_correction_single() {
        let seq = vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)];
        let num = 100;
        let mut reads: Vec<_> = (0..num).map(|_| seq.clone()).collect();
        // Introduce errors
        reads[4].insert(2, (1, 0));
        reads[5].remove(3);
        reads[2][1] = (1, 1);
        let reads: Vec<_> = reads
            .into_iter()
            .enumerate()
            .map(|(idx, seq)| (idx as u64, seq))
            .collect();
        let result = local_correction_inner(reads, (1, -1, -2), 0);
        for res in result {
            let units: Vec<_> = res.nodes.iter().map(|u| (u.unit, u.cluster)).collect();
            assert_eq!(units, seq, "{}", res.id);
        }
    }
    #[test]
    fn error_correction_multi() {
        let mut rng: Xoshiro256Plus = rand::SeedableRng::seed_from_u64(100);
        let conf = TestConfig {
            cl: 2,
            num: 100,
            fail: 0.01,
            skip: 0.01,
            max_len: 10,
            min_len: 8,
            unit_len: 40,
        };
        let reads = gen_dataset(&mut rng, conf);
        for read in reads.iter() {
            let prev = read.1.len();
            eprintln!("Correcting:{:?}", read);
            let res = correct(read, &reads, (1, -1, -2), 0);
            let seq: Vec<_> = res.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
            let after = seq.len();
            eprintln!("Corrected:{:?}", seq);
            assert!(prev.max(4) - 4 < after || after < prev + 4);
            let cl = res.nodes[0].cluster;
            let cluster = res.nodes.iter().all(|n| n.cluster == cl);
            assert!(cluster, "{:?}", res);
            let suc = res
                .nodes
                .windows(2)
                .all(|w| w[0].unit + 1 == w[1].unit || w[1].unit + 1 == w[0].unit);
            assert!(suc, "{:?}", res);
        }
    }
    #[test]
    fn unit_alignment() {
        let refr = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let query = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Match(0, 0),
                Cigar::Match(1, 0),
                Cigar::Match(2, 0),
                Cigar::Match(3, 0)
            ]
        );
        let refr = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let query = vec![(0, 0), (1, 1), (2, 0), (3, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Match(0, 0),
                Cigar::Match(1, 1),
                Cigar::Match(2, 0),
                Cigar::Match(3, 0)
            ]
        );
        let refr = vec![(1, 0), (2, 0), (3, 0)];
        let query = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Ins(0, 0),
                Cigar::Match(1, 0),
                Cigar::Match(2, 0),
                Cigar::Match(3, 0)
            ]
        );
        let refr = vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)];
        let query = vec![(0, 0), (2, 0), (3, 0), (4, 0), (5, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Match(0, 0),
                Cigar::Del,
                Cigar::Match(2, 0),
                Cigar::Match(3, 0),
                Cigar::Match(4, 0),
                Cigar::Match(5, 0)
            ]
        );
        let refr = vec![(0, 0), (1, 0), (3, 0), (4, 0), (5, 0)];
        let query = vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Match(0, 0),
                Cigar::Match(1, 0),
                Cigar::Ins(2, 0),
                Cigar::Match(3, 0),
                Cigar::Match(4, 0),
                Cigar::Match(5, 0)
            ]
        );
        let refr = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let query = vec![(0, 0), (1, 0), (2, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Match(0, 0),
                Cigar::Match(1, 0),
                Cigar::Match(2, 0),
                Cigar::Del
            ]
        );
        let refr = vec![(0, 0), (1, 0), (2, 0), (3, 0)];
        let query = vec![(1, 0), (2, 0), (3, 0), (4, 0)];
        let (_, cigar) = alignment(&query, &refr, (1, -1, -2)).unwrap();
        assert_eq!(
            cigar,
            vec![
                Cigar::Del,
                Cigar::Match(1, 0),
                Cigar::Match(2, 0),
                Cigar::Match(3, 0),
                Cigar::Ins(4, 0),
            ]
        );
    }
}

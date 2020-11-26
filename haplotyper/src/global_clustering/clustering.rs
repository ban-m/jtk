use super::error_correction::CorrectedRead;
use std::collections::HashMap;
const K: u32 = 2;
pub fn clustering_by_exact_match(reads: &[CorrectedRead]) -> (HashMap<u64, usize>, usize) {
    let graph: Vec<Vec<_>> = reads
        .iter()
        .map(|read| {
            reads
                .iter()
                .enumerate()
                .map(|(idx, target)| (idx, alignment(read, target)))
                .filter(|x| x.1 > K)
                .collect()
        })
        .collect();
    // DFS for clustering.
    let mut arrived = vec![false; graph.len()];
    let mut cluster = vec![0; graph.len()];
    let mut current_cluster = 0;
    for i in 0..graph.len() {
        if arrived[i] {
            continue;
        }
        let mut stack = vec![i];
        'dfs: while !stack.is_empty() {
            let last = *stack.last().unwrap();
            if !arrived[last] {
                arrived[last] = true;
            }
            cluster[last] = current_cluster;
            for &(to, _) in graph[last].iter() {
                if !arrived[to] {
                    stack.push(to);
                    continue 'dfs;
                }
            }
            stack.pop();
        }
        current_cluster += 1;
    }
    let mut cluster_num = vec![0; current_cluster + 1];
    for &c in cluster.iter() {
        cluster_num[c] += 1;
    }
    use std::collections::HashSet;
    let mut units: Vec<HashSet<(u64, u64)>> = vec![HashSet::new(); current_cluster + 1];
    for (i, &c) in cluster.iter().enumerate() {
        if cluster_num[c] > 10 {
            for n in reads[i].nodes.iter() {
                units[c].insert((n.unit, n.cluster));
            }
        }
    }
    let clustering: HashMap<_, _> = reads
        .iter()
        .zip(cluster)
        .filter_map(|(read, c)| {
            if cluster_num[c] > 10 {
                Some((read.id, c))
            } else {
                units
                    .iter()
                    .map(|set| {
                        read.nodes
                            .iter()
                            .filter(|n| set.contains(&(n.unit, n.cluster)))
                            .count()
                    })
                    .enumerate()
                    .max_by_key(|x| x.1)
                    .map(|(c, _)| (read.id, c))
            }
        })
        .collect();
    (clustering, current_cluster + 1)
}

fn alignment(r1: &CorrectedRead, r2: &CorrectedRead) -> u32 {
    let r1: Vec<_> = r1.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
    let r2: Vec<_> = r2.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
    let is_disjoint = {
        use std::collections::HashSet;
        let r1: HashSet<_> = r1.iter().collect();
        !r2.iter().any(|x| r1.contains(x))
    };
    if r1.is_empty() || r2.is_empty() || is_disjoint {
        return 0;
    }
    let forward = overlap_forward(&r1, &r2);
    let reverse = overlap_reverse(&r2, &r1);
    forward.max(reverse)
}

fn overlap_forward(r1: &[(u64, u64)], r2: &[(u64, u64)]) -> u32 {
    overlap_forward_inner(r1, r2).max(overlap_forward_inner(r2, r1))
}
fn overlap_reverse(r1: &[(u64, u64)], r2: &[(u64, u64)]) -> u32 {
    let (l, s) = if r1.len() > r2.len() {
        (r1, r2)
    } else {
        (r2, r1)
    };
    let longer = (0..l.len())
        .map(|x| {
            l[x.max(s.len()) - s.len()..x]
                .iter()
                .rev()
                .zip(s.iter())
                .fold((0, true), update)
                .0
        })
        .max()
        .unwrap();
    let shorter = (0..s.len())
        .map(|x| {
            s.iter()
                .rev()
                .take(x)
                .zip(l.iter().skip(l.len() - x))
                .fold((0, true), update)
                .0
        })
        .max()
        .unwrap();
    longer.max(shorter)
}

fn update((acc, exact_match): (u32, bool), (x, y): (&(u64, u64), &(u64, u64))) -> (u32, bool) {
    if x == y && exact_match {
        (acc + 1, true)
    } else {
        (0, false)
    }
}

fn overlap_forward_inner(r1: &[(u64, u64)], r2: &[(u64, u64)]) -> u32 {
    (0..r1.len())
        .map(|x| r1.iter().skip(x).zip(r2.iter()).fold((0, true), update).0)
        .max()
        .unwrap()
}

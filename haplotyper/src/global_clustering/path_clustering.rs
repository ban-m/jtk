use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
pub struct Graph {
    nodes: Vec<GNode>,
}

impl Graph {
    fn new<R: Rng>(
        reads: &[Vec<usize>],
        cs: &[usize],
        cluster: usize,
        _r: &mut R,
        remove: usize,
    ) -> Self {
        let max = *reads
            .iter()
            .flat_map(|read| read.iter().max())
            .max()
            .unwrap();
        let nodes: Vec<_> = (0..=max).map(|_| GNode::new()).collect();
        let graph = Graph { nodes };
        reads
            .iter()
            .zip(cs)
            .enumerate()
            .filter(|&(idx, (_, &cl))| cl == cluster && idx != remove)
            .fold(graph, |g, (_, (r, _))| g.register(r))
            .finalize()
    }
    fn register(mut self, read: &[usize]) -> Self {
        for w in read.windows(2) {
            self.nodes[w[0]].register(w[1]);
            self.nodes[w[1]].register(w[0]);
        }
        self
    }
    fn finalize(mut self) -> Self {
        self.nodes
            .iter_mut()
            .filter(|n| n.total > 0.0001)
            .for_each(|node| {
                let total = node.total;
                node.edges.iter_mut().for_each(|e| e.weight /= total);
            });
        self
    }
    fn score(&self, path: &[usize]) -> f64 {
        if path.len() == 0 {
            1.
        } else if path.len() == 1 {
            self.nodes[path[0]].total
        } else {
            path.windows(2)
                .map(|w| {
                    self.nodes[w[0]]
                        .edges
                        .iter()
                        .find(|e| e.to == w[1])
                        .map(|e| e.weight)
                        .unwrap_or(0.)
                })
                .sum::<f64>()
        }
    }
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let node_num = self.nodes.len();
        let edge_num = self.nodes.iter().map(|n| n.edges.len()).count();
        write!(f, "NodeNum:{}\tEdgeNum:{}", node_num, edge_num)
    }
}

struct GNode {
    edges: Vec<GEdge>,
    total: f64,
}

impl GNode {
    fn new() -> Self {
        Self {
            edges: vec![],
            total: 0.,
        }
    }
    fn register(&mut self, to: usize) {
        self.total += 1.;
        if let Some(edge) = self.edges.iter_mut().find(|e| e.to == to) {
            edge.weight += 1.;
        } else {
            self.edges.push(GEdge { to, weight: 1. });
        }
    }
}

struct GEdge {
    to: usize,
    weight: f64,
}

pub fn path_clustering(reads: &[Vec<usize>], init: &[usize]) -> Vec<usize> {
    let cluster_num = *init.iter().max().unwrap() + 1;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(100);
    let mut cluster = init.to_vec();
    assert!(!reads.is_empty());
    debug!("Start path clustering");
    let mut count = 0;
    //let stable_thr = (reads.len() / 100).max(4) as u32;
    let stable_thr = 2;
    let len = reads.len();
    while count < 3 {
        let change_num = rand::seq::index::sample(&mut rng, len, len)
            .into_iter()
            .map(|idx| {
                let models: Vec<_> = (0..cluster_num)
                    .map(|cl| Graph::new(&reads, &cluster, cl, &mut rng, idx))
                    .collect();
                let scores: Vec<_> = models.iter().map(|m| m.score(&reads[idx])).collect();
                let (argmax, _) = scores
                    .iter()
                    .enumerate()
                    .max_by(|x, y| match (x.1).partial_cmp(&y.1) {
                        Some(res) => res,
                        None => panic!("{}\t{:?}", line!(), scores),
                    })
                    .unwrap();
                let changed = 1 - (cluster[idx] == argmax) as u32;
                cluster[idx] = argmax;
                changed
            })
            .sum::<u32>();
        count += (change_num < stable_thr) as u32;
        count *= (change_num < stable_thr) as u32;
        eprintln!("ChangeNum:{},{}", change_num, stable_thr);
    }
    cluster
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
    use std::collections::{HashMap, HashSet};
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> (Vec<Vec<usize>>, Vec<usize>) {
        let TestConfig {
            cl,
            num,
            fail,
            skip,
            max_len,
            min_len,
            unit_len,
        } = conf;
        let mut answer = vec![];
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
                    let cluster = r.gen::<u64>() % (cl + 1) as u64;
                    read.push((unit as u64, cluster));
                } else {
                    read.push((unit as u64, cluster));
                }
            }
            answer.push(cluster as usize);
            reads.push(read);
        }
        let units: HashSet<_> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let units: HashMap<(u64, u64), usize> =
            units.into_iter().enumerate().map(|(x, y)| (y, x)).collect();
        let reads: Vec<Vec<_>> = reads
            .into_iter()
            .map(|read| read.into_iter().map(|x| units[&x]).collect())
            .collect();
        (reads, answer)
    }
    fn path_clustering_test() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 200,
            fail: 0.0,
            skip: 0.0,
            max_len: 20,
            min_len: 10,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let init: Vec<_> = answer
            .iter()
            .map(|&ans| {
                if ans == 0 {
                    ans
                } else if rng.gen_bool(0.5) {
                    1
                } else {
                    2
                }
            })
            .collect();
        let preds = path_clustering(&reads, &init);
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
}

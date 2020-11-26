use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let start = std::time::Instant::now();
    let ds = get_input_file()?;
    let _graph = assemble(&ds);
    let end = std::time::Instant::now();
    eprintln!("{:?}", end - start);
    Ok(())
}

fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

#[derive(Debug, Clone)]
pub struct StringGraph {
    pub nodes: Vec<StrNode>,
}

impl std::fmt::Display for StringGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let n = self.nodes.len();
        let e = self.nodes.iter().map(|e| e.edge_num()).sum::<usize>() / 2;
        write!(f, "N:{}\tE:{}", n, e)
    }
}

impl StringGraph {
    pub fn new(reads: &[&[(u64, u64)]]) -> Self {
        let nodes = reads
            .iter()
            .enumerate()
            .map(|(i, r)| StrNode::new(i, r, reads))
            .collect();
        Self { nodes }
    }
    pub fn transitive_edge_reduction(&mut self) {
        let removed_edges: Vec<Vec<bool>> = self
            .nodes
            .iter()
            .map(|node| {
                node.edges
                    .iter()
                    .enumerate()
                    .map(|(selected, e)| {
                        // Can we reach (e.to, e.to_position) from this node
                        // without using the `selected`-th edges of the node.edge?
                        let target = (e.to, e.to_position);
                        node.edges.iter().enumerate().any(|(i, f)| {
                            // If from position is different, it is not
                            // the edge we want to check.
                            if i == selected || f.from_position != e.from_position {
                                false
                            } else {
                                // Can we reach `target` node from
                                // `self.nodes[f.to]`'s `!f.to_position` node?
                                self.nodes[f.to].edges.iter().any(|g| {
                                    g.from_position == !f.to_position
                                        && (g.to, g.to_position) == target
                                })
                            }
                        })
                    })
                    .collect()
            })
            .collect();
        for (n, e) in self.nodes.iter_mut().zip(removed_edges) {
            let mut i = 0;
            n.edges.retain(|_| {
                i += 1;
                e[i - 1]
            })
        }
    }
    pub fn connected_component(&self, thr: usize) {
        use haplotyper::find_union::FindUnion;
        let mut fu = FindUnion::new(self.nodes.len());
        for node in self.nodes.iter() {
            for edge in node.edges.iter().filter(|e| e.offset > thr) {
                fu.unite(edge.to, edge.from).unwrap();
            }
        }
        let mut component = 0;
        for idx in 0..self.nodes.len() {
            if fu.find(idx).unwrap() == idx {
                let count = (0..self.nodes.len())
                    .filter(|&c| fu.find(c).unwrap() == idx)
                    .count();
                eprintln!("COMPONENT\t{}\t{}", component, count);
                component += 1;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct StrNode {
    edges: Vec<StrEdge>,
    seq: Vec<(u64, u64)>,
}

impl StrNode {
    fn edge_num(&self) -> usize {
        self.edges.len()
    }
    fn new(r_idx: usize, r: &[(u64, u64)], reads: &[&[(u64, u64)]]) -> Self {
        let mut edges = vec![];
        for (q_idx, q) in reads.iter().filter(|&&q| r != q).enumerate() {
            let indicies = (r_idx, q_idx);
            if let Some(res) = alignment((r, q), false, indicies) {
                let r: Vec<_> = r.iter().map(|x| format!("{}-{}", x.0, x.1)).collect();
                eprintln!("{}", r.join(","));
                let q: Vec<_> = q.iter().map(|x| format!("{}-{}", x.0, x.1)).collect();
                eprintln!("{}", q.join(","));
                eprintln!("{:?}\t{:?}", res.from_position, res.to_position);
                eprintln!("{}", res.offset);
                eprintln!();
                edges.push(res);
            }
            if let Some(res) = alignment((r, q), true, indicies) {
                edges.push(res)
            }
        }
        let seq = r.to_vec();
        Self { edges, seq }
    }
}

fn alignment(
    (r, q): (&[(u64, u64)], &[(u64, u64)]),
    q_forward: bool,
    indicies: (usize, usize),
) -> Option<StrEdge> {
    // Alignment
    let q = if q_forward {
        q.to_vec()
    } else {
        let mut q = q.to_vec();
        q.reverse();
        q
    };
    let r_max = (0..r.len())
        .filter_map(|r_idx| {
            if r.iter().skip(r_idx).zip(q.iter()).all(|(x, y)| x == y) {
                Some(r.len() - r_idx)
            } else {
                None
            }
        })
        .max();
    let q_max = (0..q.len())
        .filter_map(|q_idx| {
            if q.iter().skip(q_idx).zip(r.iter()).all(|(x, y)| x == y) {
                Some(q.len() - q_idx)
            } else {
                None
            }
        })
        .max();
    use Position::*;
    let (r_pos, q_pos, alnsize) = match (r_max, q_max) {
        (None, None) => return None,
        (Some(res), None) if q_forward => (Tail, Head, res),
        (Some(res), None) => (Tail, Tail, res),
        (None, Some(res)) if q_forward => (Head, Tail, res),
        (None, Some(res)) => (Head, Head, res),
        (Some(x), Some(y)) if x > y && q_forward => (Tail, Head, x),
        (Some(x), Some(y)) if x > y => (Tail, Tail, x),
        (Some(_), Some(x)) if q_forward => (Head, Tail, x),
        (Some(_), Some(x)) => (Head, Head, x),
    };
    Some(StrEdge::new(indicies, r_pos, q_pos, alnsize))
}

#[derive(Debug, Clone)]
pub struct StrEdge {
    from: usize,
    to: usize,
    from_position: Position,
    to_position: Position,
    offset: usize,
}
impl StrEdge {
    fn new(
        (r_idx, q_idx): (usize, usize),
        r_pos: Position,
        q_pos: Position,
        offset: usize,
    ) -> Self {
        Self {
            from: r_idx,
            to: q_idx,
            from_position: r_pos,
            to_position: q_pos,
            offset,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Position {
    Head,
    Tail,
}

impl std::ops::Not for Position {
    type Output = Self;
    fn not(self) -> Self::Output {
        match self {
            Self::Head => Self::Tail,
            Self::Tail => Self::Head,
        }
    }
}

impl Position {}

fn assemble(ds: &DataSet) -> StringGraph {
    let c = haplotyper::global_clustering::GlobalClusteringConfig::new(1, 1, 1, -1, -1);
    let corrected_reads = haplotyper::global_clustering::error_correction::local_correction(ds, &c);
    let reads: Vec<Vec<_>> = corrected_reads
        .iter()
        .map(|r| r.nodes.iter().map(|n| (n.unit, n.cluster)).collect())
        .collect();
    let reads: Vec<_> = reads
        .iter()
        .enumerate()
        .filter(|&(idx, ref r)| {
            let r = r.as_slice();
            reads
                .iter()
                .enumerate()
                .all(|(i, q)| i == idx || q.windows(r.len()).all(|qw| qw != r))
        })
        .map(|(_, r)| r.as_slice())
        .collect();
    let mut graph = StringGraph::new(&reads);
    eprintln!("{}", graph);
    graph.transitive_edge_reduction();
    graph
}

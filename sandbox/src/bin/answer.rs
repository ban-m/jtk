use definitions::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
#[derive(Clone, Serialize, Deserialize)]
pub struct ERead {
    pub id: u64,
    pub path: Vec<Elm>,
    pub cluster: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Elm {
    pub unit: u64,
    pub cluster: usize,
}

#[derive(Serialize, Deserialize, Clone, Copy, Hash, PartialEq, Eq)]
struct Edge {
    from: u64,
    from_cluster: usize,
    to: u64,
    to_cluster: usize,
    answer: usize,
}

#[derive(Serialize, Deserialize)]
struct Unit {
    id: usize,
    position: usize,
}
use std::fs::File;
use std::io::{BufReader, BufWriter};
// # Example
// ```bash
// cargo run --release --bin answer ${JSON(Encoded Data)} ${JSON(Log file)} > <JSON>
// ```
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let rdr = BufReader::new(File::open(&args[1])?);
    let dataset: DataSet = serde_json::de::from_reader(rdr).unwrap();
    let rdr = BufReader::new(File::open(&args[2])?);
    let ereads: Vec<ERead> = serde_json::de::from_reader(rdr).unwrap();
    let answer: Vec<usize> = {
        let max = dataset.raw_reads.iter().map(|r| r.id).max().unwrap();
        let mut answers = vec![0; max as usize + 1];
        for r in dataset.raw_reads.iter() {
            let cluster = if r.name.contains("hapA") {
                0
            } else if r.name.contains("hapB") {
                1
            } else if r.name.contains("hapC") {
                2
            } else {
                eprintln!("Unknown name:{}", r.name);
                3
            };
            answers[r.id as usize] = cluster;
        }
        answers
    };
    let (nodes, graph) = {
        use std::collections::HashSet;
        let nodes = dataset.selected_chunks.iter().map(|e| e.id).max().unwrap();
        let mut edges: Vec<_> = (0..(nodes as usize + 1)).map(|_| HashSet::new()).collect();
        for read in ereads.iter() {
            for w in read.path.windows(2) {
                let from = w[0].unit as usize;
                let to = w[1].unit as usize;
                edges[from].insert(to);
                edges[to].insert(from);
            }
        }
        let edges: Vec<Vec<_>> = edges
            .into_iter()
            .map(|eds| eds.into_iter().collect())
            .collect();
        (nodes + 1, edges)
    };
    // for (idx, eds) in graph.iter().enumerate() {
    //     eprintln!("{}\t{:?}", idx, eds);
    // }
    let sorted_units = topological_sort(&graph, nodes as usize);
    let sorted_units: HashMap<_, _> = sorted_units
        .into_iter()
        .enumerate()
        .map(|(ord, unit_id)| (unit_id, ord))
        .collect();
    let edges: Vec<(Edge, usize)> = (0..3)
        .flat_map(|cluster| {
            let mut edges: HashMap<_, usize> = HashMap::new();
            for r in ereads.iter().filter(|r| answer[r.id as usize] == cluster) {
                for w in r.path.windows(2) {
                    let x = w[0].unit;
                    let x_cluster = w[0].cluster;
                    let y = w[1].unit;
                    let y_cluster = w[1].cluster;
                    let edge = if x < y {
                        Edge {
                            from: x,
                            from_cluster: x_cluster,
                            to: y,
                            to_cluster: y_cluster,
                            answer: cluster,
                        }
                    } else {
                        Edge {
                            from: y,
                            from_cluster: y_cluster,
                            to: y,
                            to_cluster: y_cluster,
                            answer: cluster,
                        }
                    };
                    *edges.entry(edge).or_default() += 1;
                }
            }
            edges.into_iter().collect::<Vec<_>>()
        })
        .collect();
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    let sorted_units: Vec<_> = sorted_units
        .into_iter()
        .map(|(id, position)| Unit { id, position })
        .collect();
    serde_json::ser::to_writer_pretty(&mut wtr, &(sorted_units, edges)).unwrap();
    Ok(())
}

fn topological_sort(edges: &[Vec<usize>], nodes: usize) -> Vec<usize> {
    let mut stack = vec![];
    let mut arrived = vec![0; nodes];
    let mut result = vec![];
    // Search for the initial node. In other words, search a node with one degree.
    let terminal_nodes: Vec<_> = edges
        .iter()
        .enumerate()
        .filter_map(|(idx, eds)| if eds.len() <= 1 { Some(idx) } else { None })
        .collect();
    let terminal_nodes = if terminal_nodes.is_empty() {
        eprintln!("There are no tips.");
        (0..nodes).collect()
    } else {
        terminal_nodes
    };
    for start in terminal_nodes {
        if arrived[start] != 0 {
            continue;
        }
        stack.push(start);
        'dfs: while !stack.is_empty() {
            let node = *stack.last().unwrap();
            if arrived[node] == 0 {
                arrived[node] = 1;
            }
            for &to in edges[node].iter() {
                if arrived[to] == 0 {
                    stack.push(to);
                    continue 'dfs;
                }
            }
            let last = stack.pop().unwrap();
            arrived[last] = 2;
            result.push(last);
        }
    }
    assert!(arrived.iter().all(|&x| x == 2), "{:?}", arrived);
    result.reverse();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn topological_sort_test() {
        let nodes = 6;
        let edges = [vec![1], vec![0, 5], vec![4, 5], vec![], vec![2], vec![1, 2]];
        let result = topological_sort(&edges, nodes);
        assert_eq!(result, vec![3, 0, 1, 5, 2, 4]);
        let nodes = 7;
        let edges = [
            vec![1, 2],
            vec![0, 5],
            vec![0, 6],
            vec![4],
            vec![5, 3],
            vec![1, 4, 6],
            vec![5, 2],
        ];
        let result = topological_sort(&edges, nodes);
        assert_eq!(result, vec![3, 4, 5, 1, 0, 2, 6]);
    }
}

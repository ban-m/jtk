use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaploGraph {
    nodes: Vec<HapNode>,
    edges: Vec<HapEdge>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HapNode {
    id: u64,
    unit: u64,
    cluster: u64,
    parental: u32,
    maternal: u32,
    phase: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HapEdge {
    source: usize,
    target: usize,
    weight: f64,
}
use definitions::*;
use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|rdr| serde_json::de::from_reader(rdr).unwrap())?;
    let config = haplotyper::assemble::AssembleConfig::new(1, 1000, false, true, 5, 2f64);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let units = &ds.selected_chunks;
    let rt = ds.read_type;
    let mut graph =
        haplotyper::assemble::ditch_graph::DitchGraph::new(&reads, Some(units), rt, &config);
    let cov = ds.coverage.unwrap();
    let lens: Vec<_> = ds.raw_reads.iter().map(|r| r.seq().len()).collect();
    graph.assign_copy_number(cov, &lens);
    graph.remove_zero_copy_elements(&lens, 0.2);
    graph.assign_copy_number(cov, &lens);
    graph.remove_zero_copy_elements(&lens, 0.5);
    graph.assign_copy_number(cov, &lens);
    let copy_number: HashMap<_, _> = graph
        .nodes()
        .map(|node| (node.node, node.copy_number.unwrap()))
        .collect();
    let diplochunk = copy_number.values().filter(|&&c| 1 < c).count();
    log::debug!("{} nodes are non-haplochunk.", diplochunk);
    let id_is_hap1: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.contains("hapA")))
        .collect();
    //(unit,cluster)->(parental, maternal).
    let mut nodes: HashMap<_, (u32, u32)> = HashMap::new();
    let reads: Vec<(u64, Vec<_>)> = ds
        .encoded_reads
        .iter()
        .map(|r| {
            let nodes: Vec<_> = r
                .nodes
                .iter()
                .filter(|n| matches!(copy_number.get(&(n.unit,n.cluster)),Some(&x) if x <= 1 ))
                .map(|n| (n.unit, n.cluster))
                .collect();
            (r.id, nodes)
        })
        .collect();
    for (id, read) in reads.iter() {
        let is_hap1 = id_is_hap1[id];
        for &node in read.iter() {
            let occs = nodes.entry(node).or_default();
            if is_hap1 {
                occs.0 += 1;
            } else {
                occs.1 += 1;
            }
        }
    }

    // Serialize nodes.
    let nodes: Vec<HapNode> = nodes
        .iter()
        .enumerate()
        .map(|(id, (&(unit, cluster), &(parental, maternal)))| HapNode {
            id: id as u64,
            unit,
            cluster,
            parental,
            maternal,
            phase: 0,
        })
        .collect();
    // Compute index.
    let node_index: HashMap<(u64, u64), usize> = nodes
        .iter()
        .enumerate()
        .map(|(i, &HapNode { unit, cluster, .. })| ((unit, cluster), i))
        .collect();
    let mut edges: HashMap<(usize, usize), f64> = HashMap::new();
    for (_, read) in reads.iter() {
        for w in read.windows(2) {
            let (from, to) = (node_index[&w[0]], node_index[&w[1]]);
            *edges.entry((from, to)).or_default() += 1f64;
            *edges.entry((to, from)).or_default() += 1f64;
        }
    }
    for ((n1, c1), &i) in node_index.iter() {
        for ((n2, c2), &j) in node_index.iter() {
            if i != j && n1 == n2 && c1 != c2 {
                *edges.entry((i, j)).or_default() = -1.1f64;
            }
        }
    }
    // for (_, read) in reads.iter() {
    //     for (i, node) in read.iter().enumerate() {
    //         let from = node_index[node];
    //         for mode in read.iter().skip(i + 1) {
    //             let to = node_index[&mode];
    //             *edges.entry((from, to)).or_default() += 1f64;
    //             *edges.entry((to, from)).or_default() += 1f64;
    //         }
    //     }
    // }
    let edges: Vec<HapEdge> = edges
        .iter()
        .map(|(&(source, target), &weight)| HapEdge {
            source,
            target,
            weight,
        })
        .collect();
    // DUMP information.
    if 2 < args.len() {
        use std::io::{BufWriter, Write};
        let mut wtr = std::fs::File::create(&args[2]).map(BufWriter::new)?;
        let mut degrees = vec![0; nodes.len()];
        for edge in edges.iter() {
            degrees[edge.target] += 1;
            degrees[edge.source] += 1;
        }
        writeln!(&mut wtr, "unit\tcluster\tdegree\thap1\thap2")?;
        for (node, degree) in nodes.iter().zip(degrees.iter()) {
            writeln!(
                &mut wtr,
                "{}\t{}\t{}\t{}\t{}",
                node.unit, node.cluster, degree, node.parental, node.maternal
            )?;
        }
    }
    let graph = HaploGraph { nodes, edges };
    // spectral_clustering(&mut graph);
    println!("{}", serde_json::ser::to_string_pretty(&graph).unwrap());
    Ok(())
}

// fn spectral_clustering(HaploGraph { nodes, edges }: &mut HaploGraph) {
//     let num_nodes = nodes.len();
//     let edges: Vec<_> = edges
//         .iter()
//         .map(|e| (e.source, e.target, e.weight))
//         .collect();
//     let assignments = vec![];
//     for (&asn, node) in assignments.iter().zip(nodes.iter_mut()) {
//         node.phase = asn;
//     }
// }

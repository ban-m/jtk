use definitions::*;
use haplotyper::assemble::ditch_graph::Position;
use log::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let mut edges: HashMap<_, Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for w in read.nodes.windows(2) {
            let from_pos = if w[0].is_forward {
                Position::Tail
            } else {
                Position::Head
            };
            let from = (w[0].unit, w[0].cluster, from_pos);
            let to_pos = if w[1].is_forward {
                Position::Head
            } else {
                Position::Tail
            };
            let to = (w[1].unit, w[1].cluster, to_pos);
            edges.entry(from).or_default().push(to);
            edges.entry(to).or_default().push(from);
        }
    }
    edges.values_mut().for_each(|xs| {
        xs.sort();
        xs.dedup();
    });
    edges.retain(|_, vs| !vs.is_empty());
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    graph.remove_lightweight_edges(1);
    let coverage = ds.coverage.unwrap();
    let lens: Vec<_> = ds.raw_reads.iter().map(|read| read.seq().len()).collect();
    let (_, edge_cp) = graph.copy_number_estimation(coverage, &lens);
    for (&(fu, fc, fp), es) in edges.iter_mut() {
        es.retain(
            |&(tu, tc, tp)| match edge_cp.get(&(((fu, fc), fp), ((tu, tc), tp))) {
                Some(&cp) => 0 < cp,
                None => false,
            },
        );
    }
    println!("RESULT\tUnit\tCluster\tPos\tFocusUnit\tFocusCluster\tDist\tRatio");
    for node in ds.selected_chunks.iter() {
        for cl in 0..node.cluster_num as u64 {
            for pos in vec![Position::Head, Position::Tail] {
                if let Some((unit, cluster, dist, lk)) = enumerate_foci(node, cl, pos, &ds, &edges)
                {
                    println!(
                        "RESULT\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                        node.id, cl, pos, unit, cluster, dist, lk
                    );
                }
            }
        }
    }
    Ok(())
}

type Node = (u64, u64, Position);
// Return (cluster, unit, distance, lkratio), if there is putative focus.
fn enumerate_foci(
    node: &Unit,
    cluster: u64,
    pos: Position,
    ds: &DataSet,
    edges: &HashMap<Node, Vec<Node>>,
) -> Option<(u64, u64, usize, f64)> {
    let err_prob = 0.1;
    // dist_node[i] = nodes at the distance of i.
    let dist_node: Vec<Vec<_>> = enumerate_node_upto(node, cluster, pos, edges, 5)?;
    let mut lks = vec![];
    for (i, nodes) in dist_node
        .iter()
        .enumerate()
        .filter(|(_, ns)| !ns.is_empty())
    {
        debug!("{}\t{:?}", i, nodes);
        // Assume that the copy numeber is 1. Thus, each element in nodes would
        // happen in 1/nodes.len() probability.
        let mut occs: Vec<_> = vec![0; nodes.len()];
        for read in ds.encoded_reads.iter() {
            if read
                .nodes
                .iter()
                .any(|n| n.unit == node.id && n.cluster == cluster)
            {
                for (i, &(n, c)) in nodes.iter().enumerate() {
                    occs[i] += read
                        .nodes
                        .iter()
                        .filter(|node| node.unit == n && node.cluster == c)
                        .count();
                }
            }
        }
        let ith_ln = vec![-(nodes.len() as f64).ln(); nodes.len()];
        let null_prob: f64 = occs
            .iter()
            .enumerate()
            .map(|(i, &occ)| occ as f64 * ith_ln[i])
            .sum();
        let correct_lk = ((1f64 - err_prob).powi(2) + err_prob / nodes.len() as f64).ln();
        let error_lk = if nodes.len() > 1 {
            let error_in_alt = (nodes.len() as f64 - 1f64).recip();
            ((1f64 - err_prob) * err_prob * error_in_alt + err_prob / nodes.len() as f64).ln()
        } else {
            ((1f64 - err_prob) * err_prob + err_prob / nodes.len() as f64).ln()
        };
        let (alt_prob, alt_idx): (f64, _) = (0..nodes.len())
            .map(|k| {
                let lk: f64 = occs
                    .iter()
                    .map(|&x| x as f64)
                    .enumerate()
                    .map(|(i, occ)| occ * (if i == k { correct_lk } else { error_lk }))
                    .sum();
                (lk, k)
            })
            .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
            .unwrap();
        let lk_ratio = alt_prob - null_prob;
        let (n, c) = nodes[alt_idx];
        let len = nodes.len();
        let occ: usize = occs.iter().sum();
        debug!(
            "Test\t{}\t{}\t{}\t{:.2}\t{}\t{}",
            n,
            c,
            i + 1,
            lk_ratio,
            len,
            occ
        );
        lks.push((n, c, i + 1, lk_ratio));
    }
    lks.into_iter()
        .max_by(|x, y| (x.3).partial_cmp(&(y.3)).unwrap())
}

fn enumerate_node_upto(
    node: &Unit,
    cluster: u64,
    pos: Position,
    edges: &HashMap<Node, Vec<Node>>,
    radius: usize,
) -> Option<Vec<Vec<(u64, u64)>>> {
    let mut nodes_at = edges
        .get(&(node.id, cluster, pos))
        .map(|res| vec![res.clone()])?;
    if nodes_at[0].is_empty() {
        return None;
    }
    for i in 0..radius {
        let mut next_nodes = vec![];
        for &(n, c, p) in nodes_at[i].iter() {
            let moved = (n, c, !p);
            if let Some(next) = edges.get(&moved) {
                next_nodes.extend(next.iter().copied());
            }
        }
        next_nodes.sort();
        next_nodes.dedup();
        nodes_at.push(next_nodes);
    }
    let nodes_at: Vec<_> = nodes_at
        .into_iter()
        .map(|nodes| nodes.into_iter().map(|(n, c, _)| (n, c)).collect())
        .collect();
    Some(nodes_at)
}

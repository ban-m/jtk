// use sandbox::IS_MOCK;
use min_flow::convex::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let gfa: gfa::GFA = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(gfa::GFA::from_reader)?;
    let mut graph = CopyNumConvexFlowGraph::new();
    let coverage: f64 = args[2].parse().unwrap();
    let nodes: HashMap<_, _> = gfa
        .iter()
        .filter_map(|record| -> Option<_> {
            match &record.content {
                gfa::Content::Seg(x) => Some((x, &record.tags)),
                _ => None,
            }
        })
        .map(|(seg, tags)| {
            let cov: usize = tags.find("cv").unwrap().1.parse().unwrap();
            let target = cov as f64 / coverage;
            let flab = format!("{}-F", seg.sid);
            let tlab = format!("{}-T", seg.sid);
            let f_node = graph.add_node(flab);
            let t_node = graph.add_node(tlab);
            let _ = graph.add_edge(f_node, t_node, CopyFlowEdge::new(80, 0, target));
            let _ = graph.add_edge(t_node, f_node, CopyFlowEdge::new(80, 0, target));
            (seg.sid.clone(), (f_node, t_node))
        })
        .collect();
    let edges = gfa.iter().filter_map(|record| match &record.content {
        gfa::Content::Edge(edge) => Some((edge, &record.tags)),
        _ => None,
    });
    let mut is_connected: HashMap<_, _> =
        nodes.keys().map(|k| (k.clone(), (false, false))).collect();
    for (edge, tags) in edges {
        // Locate.
        if edge.beg1.is_last {
            is_connected.get_mut(&edge.sid1.id).unwrap().1 |= true;
        } else {
            is_connected.get_mut(&edge.sid1.id).unwrap().0 |= true;
        };
        if edge.beg2.is_last {
            is_connected.get_mut(&edge.sid2.id).unwrap().1 |= true;
        } else {
            is_connected.get_mut(&edge.sid2.id).unwrap().0 |= true;
        }
        let from = match edge.beg1.is_last {
            true => &nodes[&edge.sid1.id].1,
            false => &nodes[&edge.sid1.id].0,
        };
        let to = match edge.beg2.is_last {
            true => &nodes[&edge.sid2.id].1,
            false => &nodes[&edge.sid2.id].0,
        };
        let cov: usize = tags.find("cv").unwrap().1.parse().unwrap();
        let target = cov as f64 / coverage;
        graph.add_edge(*from, *to, CopyFlowEdge::new(80, 0, target));
        graph.add_edge(*to, *from, CopyFlowEdge::new(80, 0, target));
    }
    let root = graph.add_node(format!("ROOT"));
    for (key, (head_coned, tail_coned)) in is_connected {
        if !head_coned {
            let to = nodes[&key].0;
            graph.add_edge(root, to, CopyFlowEdge::mock());
            graph.add_edge(to, root, CopyFlowEdge::mock());
        }
        if !tail_coned {
            let to = nodes[&key].1;
            graph.add_edge(root, to, CopyFlowEdge::mock());
            graph.add_edge(to, root, CopyFlowEdge::mock());
        }
    }
    let flow = min_flow::min_cost_flow_convex(&graph).unwrap();
    for (eidx, cp) in flow.iter() {
        graph.edge_weight_mut(eidx).unwrap().copy = cp;
    }
    // println!("{flow}");
    use petgraph::dot::Dot;
    let edge_label = |_, _| String::new();
    let node_label = |_, _| String::new();
    let dstring = Dot::with_attr_getters(&graph, &[], &edge_label, &node_label);
    println!("{dstring:?}");
    Ok(())
}

use petgraph::graph::DiGraph;
type CopyNumConvexFlowGraph = DiGraph<String, CopyFlowEdge>;

#[derive(Clone, Copy)]
struct CopyFlowEdge {
    demand: usize,
    capacity: usize,
    target: f64,
    is_mock: bool,
    copy: usize,
}

impl std::fmt::Debug for CopyFlowEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.3}-{}", self.target, self.copy)
    }
}

impl CopyFlowEdge {
    fn new(capacity: usize, demand: usize, copy_number: f64) -> Self {
        Self {
            capacity,
            demand,
            target: copy_number,
            is_mock: false,
            copy: 0,
        }
    }
    fn mock() -> Self {
        Self {
            capacity: 80,
            demand: 0,
            target: 0f64,
            is_mock: true,
            copy: 0,
        }
    }
}

impl min_flow::flow::FlowEdge for CopyFlowEdge {
    fn demand(&self) -> min_flow::FlowRate {
        self.demand
    }
    fn capacity(&self) -> min_flow::FlowRate {
        self.capacity
    }
}

impl ConvexCost for CopyFlowEdge {
    fn convex_cost(&self, flow: min_flow::FlowRate) -> min_flow::Cost {
        match self.is_mock {
            false => (self.target - flow as f64).powi(2),
            true => 0f64,
        }
    }
}

// fn copy_number_esimation(ds:&DataSet)->HashMap<_,_> {

// }

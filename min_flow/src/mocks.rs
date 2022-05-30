use super::flow::{Flow, FlowEdgeBase, FlowGraph};
// use super::{Cost, FlowRate};
use super::FlowRate;
use petgraph::graph::{EdgeIndex, Graph};

///
/// very small flow network
/// with triangle circular topology
///
pub fn mock_flow_network1() -> (FlowGraph, Flow) {
    let mut graph: FlowGraph = Graph::new();
    let a = graph.add_node(());
    let b = graph.add_node(());
    let c = graph.add_node(());
    graph.add_edge(a, b, FlowEdgeBase::new(0, 10, -1.0));
    graph.add_edge(b, c, FlowEdgeBase::new(0, 10, -2.0));
    graph.add_edge(c, a, FlowEdgeBase::new(0, 10, -2.0));
    let f = Flow::from_vec(
        3,
        0,
        &[
            (EdgeIndex::new(0), 10),
            (EdgeIndex::new(1), 10),
            (EdgeIndex::new(2), 10),
        ],
    );
    (graph, f)
}

///
/// mock network cited from Genome-scale algorithm design p48
///
pub fn mock_flow_network2() -> (FlowGraph, Flow) {
    let mut g: FlowGraph = Graph::new();
    let s = g.add_node(());
    let a = g.add_node(());
    let b = g.add_node(());
    let c = g.add_node(());
    let d = g.add_node(());
    let e = g.add_node(());
    let f = g.add_node(());
    let t = g.add_node(());
    const INF: FlowRate = 100000;
    let e0 = g.add_edge(s, a, FlowEdgeBase::new(0, INF, 0.0));
    let e1 = g.add_edge(a, b, FlowEdgeBase::new(2, 4, 2.0));
    let e2 = g.add_edge(a, d, FlowEdgeBase::new(9, 13, 1.0));
    let e3 = g.add_edge(b, c, FlowEdgeBase::new(2, 3, 2.0));
    let e4 = g.add_edge(d, c, FlowEdgeBase::new(0, 6, 1.0));
    let e5 = g.add_edge(c, t, FlowEdgeBase::new(4, 8, 1.0));
    let e6 = g.add_edge(d, f, FlowEdgeBase::new(0, 10, 3.0));
    let e7 = g.add_edge(s, e, FlowEdgeBase::new(0, 6, 1.0));
    let e8 = g.add_edge(e, f, FlowEdgeBase::new(0, 5, 1.0));
    let e9 = g.add_edge(f, t, FlowEdgeBase::new(7, 13, 3.0));
    let ea = g.add_edge(t, s, FlowEdgeBase::new(17, 17, 0.0));

    let f = Flow::from_vec(
        11,
        0,
        &[
            (e0, 12),
            (e1, 2),
            (e2, 10),
            (e3, 2),
            (e4, 6),
            (e5, 8),
            (e6, 4),
            (e7, 5),
            (e8, 5),
            (e9, 9),
            (ea, 17),
        ],
    );

    (g, f)
}

pub fn mock_flow_network3() -> (FlowGraph, Flow) {
    let mut g: FlowGraph = Graph::new();
    let a = g.add_node(());
    let b = g.add_node(());
    let c = g.add_node(());
    let d = g.add_node(());
    let e0 = g.add_edge(a, b, FlowEdgeBase::new(0, 2, 1.0));
    let e1 = g.add_edge(a, c, FlowEdgeBase::new(0, 2, -2.0));
    let e2 = g.add_edge(b, d, FlowEdgeBase::new(0, 2, 3.0));
    let e3 = g.add_edge(c, d, FlowEdgeBase::new(0, 1, 4.0));
    let e4 = g.add_edge(d, a, FlowEdgeBase::new(2, 2, 0.0));

    let f = Flow::from_vec(5, 0, &[(e0, 1), (e1, 1), (e2, 1), (e3, 1), (e4, 2)]);
    (g, f)
}

pub fn mock_flow_network_parallel_edge1() -> (FlowGraph, Flow) {
    let mut graph: FlowGraph = Graph::new();
    let a = graph.add_node(());
    let b = graph.add_node(());
    let c = graph.add_node(());
    let e0 = graph.add_edge(a, b, FlowEdgeBase::new(0, 2, -1.0));
    let e1 = graph.add_edge(b, c, FlowEdgeBase::new(0, 2, -1.0));
    let e2 = graph.add_edge(b, c, FlowEdgeBase::new(0, 2, -2.0));
    let e3 = graph.add_edge(c, a, FlowEdgeBase::new(0, 2, 0.0));

    let f = Flow::from_vec(4, 0, &[(e0, 2), (e1, 0), (e2, 2), (e3, 2)]);
    (graph, f)
}

pub fn mock_flow_network_parallel_edge2() -> (FlowGraph, Flow) {
    let mut graph: FlowGraph = FlowGraph::new();
    let a = graph.add_node(());
    let b = graph.add_node(());
    let c = graph.add_node(());
    let e1 = graph.add_edge(a, b, FlowEdgeBase::new(0, 2, 1.0));
    let e2 = graph.add_edge(b, c, FlowEdgeBase::new(0, 2, 1.0));
    let e3 = graph.add_edge(b, c, FlowEdgeBase::new(0, 2, 2.0));
    let e4 = graph.add_edge(c, a, FlowEdgeBase::new(2, 2, 0.0));

    let f = Flow::from_vec(4, 0, &[(e1, 2), (e2, 2), (e3, 0), (e4, 2)]);
    (graph, f)
}

// #[cfg(test)]
// mod tests {
//     use super::super::utils::{draw, draw_with_flow};
//     use super::*;
//     use crate::min_flow::min_cost_flow;

//     #[test]
//     fn test_mock_flow_network1() {
//         let (g, f_true) = mock_flow_network1();
//         let f = min_cost_flow(&g).unwrap();
//         assert!(f_true == f);
//     }

//     #[test]
//     fn test_mock_flow_network2() {
//         let (g, f_true) = mock_flow_network2();
//         let f = min_cost_flow(&g).unwrap();
//         draw_with_flow(&g, &f);
//         println!("{:?}", f);
//         assert!(f_true == f);
//     }

//     #[test]
//     fn test_mock_flow_network3() {
//         let (g, f_true) = mock_flow_network3();
//         let f = min_cost_flow(&g).unwrap();
//         draw_with_flow(&g, &f);
//         println!("{:?}", f);
//         assert!(f_true == f);
//     }

//     #[test]
//     fn test_mock_flow_network_parallel_edge1() {
//         let (g, f_true) = mock_flow_network_parallel_edge1();
//         let f = min_cost_flow(&g).unwrap();
//         assert!(f_true == f);
//     }

//     #[test]
//     fn test_mock_flow_network_parallel_edge2() {
//         let (g, f_true) = mock_flow_network_parallel_edge2();
//         draw(&g);
//         let f = min_cost_flow(&g).unwrap();
//         draw_with_flow(&g, &f);

//         println!("{:?}", f);
//         println!("{:?}", f_true);
//         assert!(f_true == f);
//     }
// }

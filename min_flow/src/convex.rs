#![allow(dead_code)]
//!
//! Flow network definitions for convex cost.
//!
pub mod fast;
use super::flow::{EdgeCost, Flow, FlowEdge, FlowEdgeRaw};
use super::utils::{clamped_log, is_convex};
use super::{Cost, FlowRate};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

/// Edge of FlowGraph with convex function cost
///
/// * `cost(f)`: cost per unit flow `c(e, f)`
///
pub trait ConvexCost: FlowEdge {
    ///
    /// cost function
    /// it is a convex function of the current flow
    fn convex_cost(&self, flow: FlowRate) -> Cost;
    ///
    /// check if this edge has a finite capacity (`capacity < 100`).
    fn is_finite_capacity(&self) -> bool {
        self.capacity() < 100
    }
    ///
    /// check if the cost function of the edge is actually convex function
    fn is_convex(&self) -> bool {
        is_convex(
            |flow| self.convex_cost(flow),
            self.demand(),
            self.capacity(),
        )
    }
}

impl<E: ConvexCost> EdgeCost for E {
    fn cost(&self, flow: FlowRate) -> Cost {
        self.convex_cost(flow)
    }
}

//
// Base implementations
//

///
/// Edge attributes of ConvexFlowGraph
/// For each edge,
///   demand
///   capacity
///   convex cost function
/// are assigned.
///
#[derive(Debug, Copy, Clone)]
pub struct BaseConvexFlowEdge {
    /// demand (lower limit of flow) of the edge l(e)
    pub demand: FlowRate,
    /// capacity (upper limit of flow) of the edge u(e)
    pub capacity: FlowRate,
    /// cost function
    /// it is a convex function of the current flow
    pub convex_cost_fn: fn(FlowRate) -> Cost,
}

impl FlowEdge for BaseConvexFlowEdge {
    fn demand(&self) -> FlowRate {
        self.demand
    }
    fn capacity(&self) -> FlowRate {
        self.capacity
    }
}

impl ConvexCost for BaseConvexFlowEdge {
    fn convex_cost(&self, flow: FlowRate) -> Cost {
        (self.convex_cost_fn)(flow)
    }
}

impl BaseConvexFlowEdge {
    pub fn new(
        demand: FlowRate,
        capacity: FlowRate,
        convex_cost_fn: fn(FlowRate) -> Cost,
    ) -> BaseConvexFlowEdge {
        BaseConvexFlowEdge {
            demand,
            capacity,
            convex_cost_fn,
        }
    }
}

/// short version of BaseConvexFlowEdge::new
fn cfe(
    demand: FlowRate,
    capacity: FlowRate,
    convex_cost: fn(FlowRate) -> Cost,
) -> BaseConvexFlowEdge {
    BaseConvexFlowEdge::new(demand, capacity, convex_cost)
}

pub type ConvexFlowGraph = DiGraph<(), BaseConvexFlowEdge>;

///
///
///
#[derive(Debug, Copy, Clone)]
pub struct ConvexFlowEdgeInfo {
    /// This represents that the edge is created from the origin
    /// in ConvexFlowGraph
    origin: EdgeIndex,
    /// The FixedCostFlowEdge has demand=0 and capacity=1.
    /// if this edge has flow=1, the original edge in ConvexFlowGraph
    /// should have flow=flow_offset+1.
    flow_offset: FlowRate,
}

pub type FixedCostFlowEdge = FlowEdgeRaw<ConvexFlowEdgeInfo>;

impl FixedCostFlowEdge {
    pub fn new(
        demand: FlowRate,
        capacity: FlowRate,
        cost: Cost,
        origin: EdgeIndex,
        flow_offset: FlowRate,
    ) -> FixedCostFlowEdge {
        FixedCostFlowEdge {
            demand,
            capacity,
            cost,
            info: ConvexFlowEdgeInfo {
                origin,
                flow_offset,
            },
        }
    }
}

///
/// A special kind of flow graph, that will be converted from ConvexFlowGraph
/// Each edge has additional information of ConvexFlowEdgeInfo, that will be
/// used when convert it back to flow in the original ConvexFlowGraph.
///
pub type FixedCostFlowGraph = DiGraph<(), FixedCostFlowEdge>;

//
// conversion functions
//

///
/// convert convex flow graph to the (normal, fixed constant cost) flow graph
///
/// create an parallel edges
/// with static cost of `e.cost(f + 1) - e.cost(f)`
///
pub fn to_fixed_flow_graph<N, E>(graph: &DiGraph<N, E>) -> Option<FixedCostFlowGraph>
where
    E: FlowEdge + ConvexCost,
{
    let mut g: FixedCostFlowGraph = FixedCostFlowGraph::new();

    for e in graph.edge_indices() {
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();

        // assert that
        // (1) flow is actually convex
        if !ew.is_finite_capacity() {
            return None;
        }
        // (2) capacity is finite
        if !ew.is_convex() {
            return None;
        }

        // convert to the FixedFlowGraph
        // (1) if demand d > 0, add an edge of [demand,capacity]=[d,d]
        if ew.demand() > 0 {
            let demand = ew.demand();
            g.extend_with_edges(&[(v, w, FixedCostFlowEdge::new(demand, demand, 0.0, e, 0))]);
        }
        // (2) add aux edges of [demand,capacity]=[0,1]
        let edges: Vec<(NodeIndex, NodeIndex, FixedCostFlowEdge)> = (ew.demand()..ew.capacity())
            .map(|f| {
                let cost = ew.cost(f + 1) - ew.cost(f);
                let fe = FixedCostFlowEdge::new(0, 1, cost, e, f);
                (v, w, fe)
            })
            .collect();
        g.extend_with_edges(&edges);
    }

    Some(g)
}

fn get_flow_in_fixed(
    edge: EdgeIndex,
    fixed_flow: &Flow,
    fixed_graph: &FixedCostFlowGraph,
) -> FlowRate {
    // for original edge e (with EdgeIndex edge) in ConvexFlowGraph
    // the flow is the sum of the flow on the edges whose FixedCostFlowEdge.info.origin == edge
    fixed_graph
        .edge_indices()
        .filter_map(|fe| {
            let fe_weight = fixed_graph.edge_weight(fe).unwrap();
            if fe_weight.info.origin == edge {
                let flow = fixed_flow[fe];
                Some(flow)
            } else {
                None
            }
        })
        .sum()
}

pub fn restore_convex_flow<N, E>(
    fixed_flow: &Flow,
    fixed_graph: &FixedCostFlowGraph,
    graph: &DiGraph<N, E>,
) -> Flow
where
    E: FlowEdge + ConvexCost,
{
    let mut flow = Flow::new(graph.edge_count(), 0);

    for e in graph.edge_indices() {
        let f = get_flow_in_fixed(e, fixed_flow, fixed_graph);
        flow[e] = f;
    }

    flow
}

//
// utils
//
#[allow(dead_code)]
fn mock_convex_flow_graph1() -> (ConvexFlowGraph, Flow) {
    let mut g: ConvexFlowGraph = ConvexFlowGraph::new();
    let a = g.add_node(());
    let b = g.add_node(());
    let c = g.add_node(());
    let e0 = g.add_edge(a, b, cfe(0, 10, |f| (f as f64 - 5.0).powi(2)));
    let e1 = g.add_edge(b, c, cfe(0, 10, |f| (f as f64 - 5.0).powi(2)));
    let e2 = g.add_edge(c, a, cfe(0, 10, |f| (f as f64 - 5.0).powi(2)));

    let f = Flow::from_vec(3, 0, &[(e0, 5), (e1, 5), (e2, 5)]);

    (g, f)
}

#[allow(dead_code)]
fn mock_convex_flow_graph2() -> (ConvexFlowGraph, Flow) {
    let mut g: ConvexFlowGraph = ConvexFlowGraph::new();
    let s = g.add_node(());
    let v1 = g.add_node(());
    let v2 = g.add_node(());
    let w1 = g.add_node(());
    let w2 = g.add_node(());
    let t = g.add_node(());
    // surrounding
    let e0 = g.add_edge(s, v1, cfe(2, 2, |_| 0.0));
    let e1 = g.add_edge(s, v2, cfe(4, 4, |_| 0.0));
    let e2 = g.add_edge(w1, t, cfe(1, 1, |_| 0.0));
    let e3 = g.add_edge(w2, t, cfe(5, 5, |_| 0.0));
    let e4 = g.add_edge(t, s, cfe(6, 6, |_| 0.0));

    // intersecting
    let e5 = g.add_edge(v1, w1, cfe(0, 6, |f| -10.0 * clamped_log(f)));
    let e6 = g.add_edge(v1, w2, cfe(0, 6, |f| -10.0 * clamped_log(f)));
    let e7 = g.add_edge(v2, w1, cfe(0, 6, |f| -10.0 * clamped_log(f)));
    let e8 = g.add_edge(v2, w2, cfe(0, 6, |f| -10.0 * clamped_log(f)));

    // true flow
    let f = Flow::from_vec(
        9,
        0,
        &[
            (e0, 2),
            (e1, 4),
            (e2, 1),
            (e3, 5),
            (e4, 6),
            (e5, 0),
            (e6, 2),
            (e7, 1),
            (e8, 3),
        ],
    );

    (g, f)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::utils::{draw, draw_with_flow};
    use super::*;
    use crate::{min_cost_flow, min_cost_flow_convex, min_cost_flow_convex_fast};

    #[test]
    fn convex_flow_edge_new() {
        let e = BaseConvexFlowEdge::new(1, 5, |f| 10.0 * f as f64);
        assert_eq!(50.0, e.convex_cost(5));
        assert_eq!(50.0, e.cost(5));
        assert_eq!(1, e.demand);
        assert_eq!(5, e.capacity);
        assert_eq!(true, e.is_convex());
    }

    #[test]
    fn convex_flow_graph_mock1() {
        let (g, f_true) = mock_convex_flow_graph1();
        let fg = to_fixed_flow_graph(&g).unwrap();
        let fg_flow = min_cost_flow(&fg).unwrap();
        let flow = restore_convex_flow(&fg_flow, &fg, &g);

        println!("{:?}", flow);
        println!("{:?}", f_true);
        assert!(flow == f_true);
    }

    #[test]
    fn convex_flow_graph_mock2() {
        let (g, f_true) = mock_convex_flow_graph2();
        draw(&g);
        let flow = min_cost_flow_convex(&g).unwrap();
        draw_with_flow(&g, &flow);
        println!("{:?}", flow);
        assert!(flow == f_true);
    }

    #[test]
    fn convex_flow_graph_mock1_fast() {
        let (g, f_true) = mock_convex_flow_graph1();
        let f = min_cost_flow_convex_fast(&g);
        println!("{:?}", f);
        println!("{:?}", f_true);
        assert!(f.is_some());
        assert!(f.unwrap() == f_true);
    }

    #[test]
    fn convex_flow_graph_mock2_fast() {
        let (g, f_true) = mock_convex_flow_graph2();
        let f = min_cost_flow_convex_fast(&g);
        println!("{:?}", f);
        println!("{:?}", f_true);
        assert!(f.is_some());
        assert!(f.unwrap() == f_true);
    }
}

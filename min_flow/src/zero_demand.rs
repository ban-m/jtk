//!
//! Zero demand flow graphs
//! for finding init valid flow
//!
use super::flow::{total_cost, Flow, FlowEdge, FlowEdgeRaw};
use super::min_cost_flow_from_zero;
use super::{Cost, FlowRate};
use petgraph::graph::{DiGraph, EdgeIndex};

// basic definitions

#[allow(dead_code)]
#[derive(Debug, Copy, Clone)]
pub struct ZeroDemandEdgeInfo {
    /// edge id in original graph
    origin: EdgeIndex,
    /// identifier of type-A or type-B
    kind: ZeroDemandEdgeKind,
}

#[derive(Debug, Copy, Clone)]
pub enum ZeroDemandEdgeKind {
    /// represents type-A edge
    /// e = ([0,l],-1)
    BelowDemand,
    /// represents type-B edge
    /// e = ([0,u-l],0)
    AboveDemand,
}

pub type ZeroDemandFlowEdge = FlowEdgeRaw<ZeroDemandEdgeInfo>;

impl ZeroDemandFlowEdge {
    pub fn new(
        demand: FlowRate,
        capacity: FlowRate,
        cost: Cost,
        origin: EdgeIndex,
        kind: ZeroDemandEdgeKind,
    ) -> ZeroDemandFlowEdge {
        ZeroDemandFlowEdge {
            demand,
            capacity,
            cost,
            info: ZeroDemandEdgeInfo { origin, kind },
        }
    }
}
pub type ZeroDemandFlowGraph = DiGraph<(), ZeroDemandFlowEdge>;

// conversion functions

/// Convert normal FlowGraph to ZeroDemandGraph
/// (min-cost-flow on ZeroDemandGraph) == (one of the valid flow on FlowGraph)
fn to_zero_demand_graph<N, E: FlowEdge + std::fmt::Debug>(
    graph: &DiGraph<N, E>,
) -> ZeroDemandFlowGraph {
    let mut zdg: ZeroDemandFlowGraph = ZeroDemandFlowGraph::new();

    // create two edges (type-A and type-B) for each edge
    for e in graph.edge_indices() {
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();
        // println!("{:?} {:?}", e, ew);

        zdg.extend_with_edges(&[
            (
                v,
                w,
                ZeroDemandFlowEdge::new(0, ew.demand(), -1., e, ZeroDemandEdgeKind::BelowDemand),
            ),
            (
                v,
                w,
                ZeroDemandFlowEdge::new(
                    0,
                    ew.capacity() - ew.demand(),
                    0.,
                    e,
                    ZeroDemandEdgeKind::AboveDemand,
                ),
            ),
        ]);
    }
    zdg
}

/// Convert-back to the original flow
/// Given the converted zero-demand graph and a flow on it (zd_graph, zd_flow),
/// this calculates the original flow on the original graph.
/// For each edge $e$ in the original graph,
/// there is
/// - type-A edge $ea$
/// - type-B edge $eb$
/// and the sum of the flows of the two is the flow of original edge $e$.
fn zero_demand_flow_to_original_flow<N, E: FlowEdge>(
    graph: &DiGraph<N, E>,
    zd_flow: &Flow,
    zd_graph: &ZeroDemandFlowGraph,
) -> Flow {
    let mut flow = Flow::new(graph.edge_count(), 0);
    for e in zd_graph.edge_indices() {
        let ew = zd_graph.edge_weight(e).unwrap();
        let zd_f = zd_flow[e];
        let original_edge = ew.info.origin;
        flow[original_edge] += zd_f;
    }
    flow
}

// public function

///
/// Find initial flow of the FlowGraph, by
/// 1. convert to ZeroDemandFlowGraph
/// 2. find the min-cost-flow on the zero demand graph
/// 3. convert it back to the valid flow on the original graph
///
pub fn find_initial_flow<N, E: FlowEdge + std::fmt::Debug>(graph: &DiGraph<N, E>) -> Option<Flow> {
    let zdg = to_zero_demand_graph(graph);
    // utils::draw(&zdg);
    let zd_flow = min_cost_flow_from_zero(&zdg);

    // println!("sum_of_demand={:?}", sum_of_demand(&graph));
    if total_cost(&zdg, &zd_flow) > sum_of_demand(&graph) as Cost {
        // valid flow does not exists
        None
    } else {
        // valid initial flow exists!
        let flow = zero_demand_flow_to_original_flow(graph, &zd_flow, &zdg);
        Some(flow)
    }
}

// utils

///
/// To determine all-zero flow is valid or not
/// we should know whether the given graph is demand-less
/// that is all demands of the edges are 0.
///
pub fn is_zero_demand_flow_graph<N, E: FlowEdge>(graph: &DiGraph<N, E>) -> bool {
    graph.edge_indices().all(|e| {
        let ew = graph.edge_weight(e).unwrap();
        ew.demand() == 0
    })
}

///
/// sum of edge demand
///
fn sum_of_demand<N, E: FlowEdge>(graph: &DiGraph<N, E>) -> FlowRate {
    graph
        .edge_indices()
        .map(|e| {
            let ew = graph.edge_weight(e).unwrap();
            ew.demand()
        })
        .sum()
}

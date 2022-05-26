#![feature(vec_retain_mut)]
// Originally written by Ryo Nakabayashi, 2022.
// Modified by Bansho Masutani, 2022.

pub mod common;
pub mod convex;
pub mod flow;
pub mod mocks;
pub mod prob;
pub mod residue;
pub mod utils;
pub mod vector;
pub mod zero_demand;
use convex::{restore_convex_flow, to_fixed_flow_graph, ConvexCost};
pub use flow::total_cost;
use flow::{is_valid_flow, ConstCost, Flow, FlowEdge};
use petgraph::graph::DiGraph;
use residue::{improve_flow, improve_flow_convex};
// use utils::draw_with_flow;
use zero_demand::{find_initial_flow, is_zero_demand_flow_graph};

/// type of a flow (on edges) in min-flow.
pub type FlowRate = usize;

/// type of a cost (of edges per unit flow) in min-flow.
pub type Cost = f64;

//
// public functions
//

///
/// Find minimum cost flow on the FlowGraph
///
pub fn min_cost_flow<N, E>(graph: &DiGraph<N, E>) -> Option<Flow>
where
    N: std::fmt::Debug,
    E: FlowEdge + ConstCost + std::fmt::Debug,
{
    let init_flow = find_initial_flow(graph);
    match init_flow {
        Some(flow) => {
            // draw_with_flow(graph, &flow);
            Some(min_cost_flow_from(graph, &flow))
        }
        None => None,
    }
}

///
/// Find minimum cost flow on the ConvexFlowGraph
///
pub fn min_cost_flow_convex<N, E>(graph: &DiGraph<N, E>) -> Option<Flow>
where
    N: std::fmt::Debug,
    E: FlowEdge + ConvexCost + std::fmt::Debug,
{
    // (1) convert to normal FlowGraph and find the min-cost-flow
    let fg = match to_fixed_flow_graph(graph) {
        Some(fg) => fg,
        None => return None,
    };
    let fg_flow = match min_cost_flow(&fg) {
        Some(fg_flow) => fg_flow,
        None => return None,
    };
    // (2) convert-back to the flow on the ConvexFlowGraph
    Some(restore_convex_flow(&fg_flow, &fg, &graph))
}

///
/// Find minimum cost flow on the Graph whose edge is ConvexFlowEdge.
/// This solver requires less memory.
///
pub fn min_cost_flow_convex_fast<N, E>(graph: &DiGraph<N, E>) -> Option<Flow>
where
    N: std::fmt::Debug,
    E: FlowEdge + ConvexCost + std::fmt::Debug,
{
    // (1) find the initial flow, by assigning constant cost to the flow.
    let init_flow = find_initial_flow(graph);

    // (2) upgrade the flow, by finding a negative cycle in residue graph.
    match init_flow {
        Some(flow) => {
            // draw_with_flow(graph, &flow);
            Some(min_cost_flow_from_convex(graph, &flow))
        }
        None => None,
    }
}

//
// internal functions
//

///
/// Find minimum cost flow of the special FlowGraph, whose demand is always zero.
///
fn min_cost_flow_from_zero<N, E: FlowEdge + ConstCost>(graph: &DiGraph<N, E>) -> Flow {
    assert!(is_zero_demand_flow_graph(&graph));
    let flow = Flow::new(graph.edge_count(), 0);
    min_cost_flow_from(graph, &flow)
}

///
/// Find minimum cost by starting from the specified flow values.
///
fn min_cost_flow_from<N, E: FlowEdge + ConstCost>(graph: &DiGraph<N, E>, init_flow: &Flow) -> Flow {
    let mut flow = init_flow.clone();

    loop {
        assert!(is_valid_flow(&flow, &graph));
        match improve_flow(graph, &flow) {
            Some(new_flow) => {
                flow = new_flow;
                continue;
            }
            None => {
                break;
            }
        };
    }

    flow
}

///
/// Find minimum cost by starting from the specified flow values in ConvexCost Flowgraph.
///
fn min_cost_flow_from_convex<N, E: FlowEdge + ConvexCost>(
    graph: &DiGraph<N, E>,
    init_flow: &Flow,
) -> Flow {
    let mut flow = init_flow.clone();

    loop {
        assert!(is_valid_flow(&flow, &graph));
        match improve_flow_convex(graph, &flow) {
            Some(new_flow) => {
                flow = new_flow;
                continue;
            }
            None => {
                break;
            }
        };
    }

    flow
}

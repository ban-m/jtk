//! Residue graph related definitions
//! - ResidueEdge
//! - ResidueGraph
//! - ResidueDirection
//!
use super::convex::ConvexCost;
use super::flow::{ConstCost, EdgeCost, Flow, FlowEdge};
// use super::utils::draw;
use super::{Cost, FlowRate};
use itertools::Itertools; // for tuple_windows
use petgraph::algo::find_negative_cycle;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::prelude::*;
use petgraph::visit::VisitMap;
use std::cmp::Ordering;

// basic definitions

/// Edge attributes used in ResidueGraph
#[derive(Debug, Default, Copy, Clone)]
pub struct ResidueEdge {
    /// The movable amount of the flow
    pub count: FlowRate,
    /// Cost of the unit change of this flow
    pub weight: Cost,
    /// Original edge index of the source graph
    pub target: EdgeIndex,
    /// +1 or -1
    pub direction: ResidueDirection,
}

impl ResidueEdge {
    pub fn new(
        count: FlowRate,
        weight: Cost,
        target: EdgeIndex,
        direction: ResidueDirection,
    ) -> ResidueEdge {
        ResidueEdge {
            count,
            weight,
            target,
            direction,
        }
    }
    pub fn only_weight(weight: Cost) -> ResidueEdge {
        ResidueEdge {
            weight,
            // filled by default values
            target: EdgeIndex::new(0),
            count: 0,
            direction: ResidueDirection::Up,
        }
    }
}

// Implement FloatMeasure for ResidueEdge
// to use ResidueEdge.weight as a weight in bellman ford
impl PartialOrd for ResidueEdge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.weight.partial_cmp(&other.weight)
    }
}

impl PartialEq for ResidueEdge {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
    }
}

impl std::ops::Add for ResidueEdge {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ResidueEdge::only_weight(self.weight + other.weight)
    }
}

impl petgraph::algo::FloatMeasure for ResidueEdge {
    fn zero() -> Self {
        ResidueEdge::only_weight(0.)
    }
    fn infinite() -> Self {
        ResidueEdge::only_weight(1. / 0.)
    }
}

/// Residue direction enum
/// residue edge has two types
#[derive(Debug, Copy, Clone)]
pub enum ResidueDirection {
    /// Up edge: it can increase(+1) flow
    Up,
    /// Down edge: it can decrease(-1) flow
    Down,
}

impl Default for ResidueDirection {
    fn default() -> Self {
        ResidueDirection::Up
    }
}

/// ResidueGraph definition
pub type ResidueGraph = DiGraph<(), ResidueEdge>;

//
// conversion functions
//

/// Convert FlowGraph with Flow into ResidueGraph.
///
/// FlowGraph and Flow
/// v -> w
///  e = ([l,u],c), f
///
/// into
///
/// ResidueGraph
/// v -> w
///  e1 = (u-f, +c) if u-f>0
/// w -> v
///  e2 = (f-l, -c) if f-l>0
pub fn flow_to_residue<N, E: FlowEdge + ConstCost>(
    graph: &DiGraph<N, E>,
    flow: &Flow,
) -> ResidueGraph {
    let mut rg: ResidueGraph = ResidueGraph::new();

    // create two edges (Up and Down) for each edge
    for e in graph.edge_indices() {
        let f = flow[e];
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();

        let mut edges = Vec::new();
        if f < ew.capacity() {
            // up movable
            edges.push((
                v,
                w,
                ResidueEdge::new(ew.capacity() - f, ew.cost(), e, ResidueDirection::Up),
            ));
        }
        if f > ew.demand() {
            // down movable
            edges.push((
                w,
                v,
                ResidueEdge::new(f - ew.demand(), -ew.cost(), e, ResidueDirection::Down),
            ));
        }
        rg.extend_with_edges(&edges);
    }
    rg
}

/// Convert FlowGraph with Flow with ConvexCost into ResidueGraph.
///
/// For each edge in FlowGraph with Flow
/// ```text
/// e(v -> w) = ([l,u],c), f
/// ```
///
/// create two edges in ResidueGraph
/// ```text
/// e1(v -> w) = (1, c(f+1) - c(f)) if u - f > 0
///
/// e2(w -> v) = (1, c(f-1) - c(f)) if f - l > 0
/// ```
pub fn flow_to_residue_convex<N, E>(graph: &DiGraph<N, E>, flow: &Flow) -> ResidueGraph
where
    E: FlowEdge + ConvexCost,
{
    let mut rg: ResidueGraph = ResidueGraph::new();

    // create two edges (Up and Down) for each edge
    for e in graph.edge_indices() {
        let f = flow[e];
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();

        let mut edges = Vec::new();
        if f < ew.capacity() {
            // up movable
            edges.push((
                v,
                w,
                ResidueEdge::new(1, ew.cost(f + 1) - ew.cost(f), e, ResidueDirection::Up),
            ));
        }
        if f > ew.demand() {
            // down movable
            edges.push((
                w,
                v,
                ResidueEdge::new(1, ew.cost(f - 1) - ew.cost(f), e, ResidueDirection::Down),
            ));
        }
        rg.extend_with_edges(&edges);
    }
    rg
}

#[allow(dead_code)]
fn residue_to_float_weighted_graph(graph: &ResidueGraph) -> DiGraph<(), Cost> {
    graph.map(|_, _| (), |_, ew| ew.weight)
}

//
// internal functions to find a update of the flow
// (i.e. the negative cycle in ResidueGraph)
//

/// Find the minimum weight edge among all parallel edges between v and w
/// Input: two nodes (v,w) in a graph
/// Output: minimum weight edge among all parallel edge (v,w)
///
/// (Used in `node_list_to_edge_list`)
fn pick_minimum_weight_edge(graph: &ResidueGraph, v: NodeIndex, w: NodeIndex) -> EdgeIndex {
    let er = graph
        .edges_connecting(v, w)
        .min_by(|e1, e2| {
            let w1 = e1.weight().weight;
            let w2 = e2.weight().weight;
            w1.partial_cmp(&w2).unwrap()
        })
        .unwrap();
    let e = er.id();
    e
}

/// Convert "a cycle as nodes [NodeIndex]" into "a cycle as edges [EdgeIndex]",
/// by choosing the minimum weight edge if there are parallel edges
fn node_list_to_edge_list(graph: &ResidueGraph, nodes: &[NodeIndex]) -> Vec<EdgeIndex> {
    let mut edges = Vec::new();

    // (1) nodes[i] and nodes[i+1] from i=0..n-1
    for (v, w) in nodes.iter().tuple_windows() {
        let edge = pick_minimum_weight_edge(graph, *v, *w);
        edges.push(edge);
    }

    // (2) tail and head, node[n-1] and node[0]
    let edge = pick_minimum_weight_edge(graph, nodes[nodes.len() - 1], nodes[0]);
    edges.push(edge);

    edges
}

fn is_negative_cycle(graph: &ResidueGraph, edges: &[EdgeIndex]) -> bool {
    let total_weight: Cost = edges
        .iter()
        .map(|&e| {
            let ew = graph.edge_weight(e).unwrap();
            ew.weight
        })
        .sum();
    total_weight < 0.0
}

///
/// Update the flow by a negative cycle on a residue graph.
///
fn apply_residual_edges_to_flow(flow: &Flow, rg: &ResidueGraph, edges: &[EdgeIndex]) -> Flow {
    let mut new_flow = flow.clone();

    // (1) determine flow_change_amount
    // that is the minimum of ResidueEdge.count
    let flow_change_amount = edges
        .iter()
        .map(|&e| {
            let ew = rg.edge_weight(e).unwrap();
            ew.count
        })
        .min()
        .unwrap();

    // (2) apply these changes to the flow along the cycle
    for edge in edges {
        let ew = rg.edge_weight(*edge).unwrap();
        // convert back to the original edgeindex
        let original_edge = ew.target;

        new_flow[original_edge] = match ew.direction {
            ResidueDirection::Up => flow[original_edge] + flow_change_amount,
            ResidueDirection::Down => flow[original_edge] - flow_change_amount,
        };
    }

    new_flow
}

fn find_negative_cycle_in_whole_graph(graph: &ResidueGraph) -> Option<Vec<NodeIndex>> {
    let mut node = NodeIndex::new(0);
    let mut dfs = Dfs::new(&graph, node);

    loop {
        let path = find_negative_cycle(&graph, node);

        if path.is_some() {
            return path;
        }

        // search for alternative start point
        dfs.move_to(node);
        while let Some(_nx) = dfs.next(&graph) {}
        let unvisited_node = graph
            .node_indices()
            .find(|node| !dfs.discovered.is_visited(node));

        // if there is unvisited node, search again for negative cycle
        match unvisited_node {
            Some(n) => {
                node = n;
                continue;
            }
            None => break,
        }
    }
    return None;
}

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
fn update_flow_in_residue_graph(flow: &Flow, rg: &ResidueGraph) -> Option<Flow> {
    // find negative weight cycles
    let path = find_negative_cycle_in_whole_graph(&rg);
    // draw(&rg);

    match path {
        Some(nodes) => {
            let edges = node_list_to_edge_list(&rg, &nodes);

            // check if this is actually negative cycle
            assert!(is_negative_cycle(&rg, &edges));

            // apply these changes along the cycle to current flow
            let new_flow = apply_residual_edges_to_flow(&flow, &rg, &edges);
            // println!("{:?}", new_flow);
            Some(new_flow)
        }
        None => None,
    }
}

//
// public functions
//

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow<N, E: FlowEdge + ConstCost>(
    graph: &DiGraph<N, E>,
    flow: &Flow,
) -> Option<Flow> {
    let rg = flow_to_residue(graph, flow);
    update_flow_in_residue_graph(flow, &rg)
}

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow_convex<N, E>(graph: &DiGraph<N, E>, flow: &Flow) -> Option<Flow>
where
    E: FlowEdge + ConvexCost,
{
    let rg = flow_to_residue_convex(graph, flow);
    update_flow_in_residue_graph(flow, &rg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn petgraph_negative_cycle_test() {
        // small cycle test
        let mut g: DiGraph<(), f32> = Graph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        g.add_edge(a, b, -10.0);
        g.add_edge(b, a, 9.0);
        let path = find_negative_cycle(&g, NodeIndex::new(0));
        assert_eq!(path.is_some(), true);
        let nodes = path.unwrap();
        assert!(nodes.contains(&NodeIndex::new(0)));
        assert!(nodes.contains(&NodeIndex::new(1)));
    }

    #[test]
    fn petgraph_negative_cycle_test2() {
        // self loop test, it will work fine
        let mut g: DiGraph<(), f32> = Graph::new();
        let a = g.add_node(());
        g.add_edge(a, a, -10.0);
        let path = find_negative_cycle(&g, NodeIndex::new(0));
        assert_eq!(path.is_some(), true);
        let nodes = path.unwrap();
        assert!(nodes.contains(&NodeIndex::new(0)));
    }

    #[test]
    fn negative_cycle_in_whole() {
        let mut g: ResidueGraph = ResidueGraph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.add_edge(
            a,
            b,
            ResidueEdge::new(1, 10.0, EdgeIndex::new(0), ResidueDirection::Up),
        );
        g.add_edge(
            b,
            a,
            ResidueEdge::new(1, -1.0, EdgeIndex::new(0), ResidueDirection::Up),
        );
        g.add_edge(
            c,
            c,
            ResidueEdge::new(1, -1.0, EdgeIndex::new(0), ResidueDirection::Up),
        );
        let path = find_negative_cycle_in_whole_graph(&g);
        assert_eq!(path.is_some(), true);
        let nodes = path.unwrap();
        assert!(nodes.contains(&NodeIndex::new(2)));
    }
}

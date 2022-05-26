//!
//! Fast min-flow solver for convex cost functions
//!
//! # (1) naive
//!
//! * convert to the constant cost graph, by duplicating edges which has a constant cost of `f(i+1)-f(i)`.
//! * solve the graph by normal min-flow solver.
//!
//! # (2) fast
//!
//! * to find the initial valid flow, cost will be unit and...
//! * to find the min flow from the init flow, residual...
//!
// use super::super::flow::{Flow, FlowEdge};
// use super::ConvexCost;
// use petgraph::graph::DiGraph;

//!
//! utils
//!
use super::flow::Flow;
use super::mocks;
use super::{find_initial_flow, min_cost_flow};
// use super::{Cost, FlowRate};
use super::FlowRate;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, Graph};
use petgraph::EdgeType;

///
/// check if the function `f` is convex or not
/// in the domain `[x_min, x_max]`
///
/// it will check `f(x + 1) - f(x)` is monotonically increasing
/// for increasing `x`s
///
pub fn is_convex<F: Fn(usize) -> f64>(f: F, x_min: usize, x_max: usize) -> bool {
    let mut y_prev = f64::MIN;

    (x_min..x_max).map(|x| f(x + 1) - f(x)).all(|y| {
        // check if ys are increasing
        let is_increasing = y >= y_prev;
        y_prev = y;
        is_increasing
    })
}

/// Avoid `ln(0) = -\infty`
/// by setting ln(0) = (some constant < 0)
pub fn clamped_log(x: usize) -> f64 {
    if x == 0 {
        -1000.0
    } else {
        (x as f64).ln()
    }
}

pub fn test() {
    let (g, _) = mocks::mock_flow_network2();
    draw(&g);

    let f = find_initial_flow(&g);
    println!("initia_flow={:?}", f);

    let f = min_cost_flow(&g);
    println!("{:?}", f);
}

pub fn draw<'a, N: 'a, E: 'a, Ty, Ix>(graph: &'a Graph<N, E, Ty, Ix>)
where
    E: std::fmt::Debug,
    N: std::fmt::Debug,
    Ty: EdgeType,
    Ix: petgraph::graph::IndexType,
{
    println!("{:?}", Dot::with_config(&graph, &[]));
}

#[allow(dead_code)]
#[derive(Debug, Copy, Clone)]
struct EdgeWithFlow<T> {
    flow: FlowRate,
    info: T,
}

pub fn draw_with_flow<N, E>(graph: &DiGraph<N, E>, flow: &Flow)
where
    N: std::fmt::Debug,
    E: std::fmt::Debug,
{
    let graph_with_flow = graph.map(
        |_, vw| vw,
        |e, ew| EdgeWithFlow {
            flow: flow[e],
            info: ew,
        },
    );
    println!("{:?}", Dot::with_config(&graph_with_flow, &[]));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_convex_test() {
        // f(x) = (x-10)^2
        assert!(is_convex(|x| (x as f64 - 10.0).powi(2), 0, 20));
        // f(x) = - (x-10)^2
        assert!(!is_convex(|x| -(x as f64 - 10.0).powi(2), 0, 20));
        // f(x) = 0 (constant)
        assert!(is_convex(|_| 0.0, 0, 20));
        // f(x) = -c log(x)
        assert!(is_convex(|x| -10.0 * (x as f64).ln(), 1, 20));
        // f(x) = -c clamped_log(x)
        assert!(is_convex(|x| -10.0 * clamped_log(x), 0, 20));
    }
}

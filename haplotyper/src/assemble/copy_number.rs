// Each edge has its direction to either "plus" or "minus" direction of the node.
type Edge = (usize, bool, usize, bool, f64);

#[derive(Debug, Clone)]
pub struct CoverageCalibrator {
    // Sorted.
    lengths: Vec<usize>,
    // i->sum of the inverse value from i to the end.
    cum_inverse_sum: Vec<f64>,
}

impl CoverageCalibrator {
    pub fn new(lens: &[usize]) -> Self {
        let mut lengths = lens.to_vec();
        lengths.sort_unstable();
        let (mut cum_inverse_sum, _): (Vec<_>, _) = lengths
            .iter()
            .map(|&x| (x as f64).recip())
            .rev()
            .fold((vec![], 0f64), |(mut sums, cum), x| {
                sums.push(cum + x);
                (sums, cum + x)
            });
        cum_inverse_sum.reverse();
        Self {
            lengths,
            cum_inverse_sum,
        }
    }
    pub fn calib(&self, observed: usize, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) => idx,
        };
        (observed * self.lengths.len()) as f64
            / ((self.lengths.len() - idx) as f64 - gap_len as f64 * self.cum_inverse_sum[idx])
    }
    pub fn calib_f64(&self, observed: f64, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) => idx,
        };
        observed * self.lengths.len() as f64
            / ((self.lengths.len() - idx) as f64 - gap_len as f64 * self.cum_inverse_sum[idx])
    }
}

// Optimizer.
#[derive(Debug, Clone)]
struct Optimizer {
    nodes: Vec<f64>,
    edges: Vec<Edge>,
    is_tip: Vec<bool>,
    gradient: Vec<f64>,
    momentum: Vec<f64>,
}

impl Optimizer {
    fn new(nodes: &[f64], edges: &[Edge], is_tip: &[bool]) -> Self {
        Self {
            nodes: nodes.to_vec(),
            edges: edges.to_vec(),
            is_tip: is_tip.to_vec(),
            gradient: vec![0f64; edges.len()],
            momentum: vec![0f64; edges.len()],
        }
    }
    fn optimize(&mut self, copy_num: &mut [f64]) {
        //debug!("LOSS\tEdgeRes\tNodeRes\tNodePenalty\tIntPenalty\tTotal");
        use std::collections::VecDeque;
        let heat = (1..=100).map(|x| (1f64, 0.01 * x as f64, 0.01 * x as f64));
        let chill = (1..=100).map(|x| (1f64 - 0.01 * x as f64, 1f64, 1f64));
        let re_heat = (1..=100).map(|x| (0.01 * x as f64, 1f64, 1f64));
        let re_chill = (1..=100).map(|x| (1f64 - 0.01 * x as f64, 1f64, 1f64 - 0.01 * x as f64));
        let finalize = (1..=50).map(|x| (0.1, 1f64 - 0.01 * x as f64, 0.01 * x as f64));
        let schedule = heat
            .chain(chill)
            .chain(re_heat)
            .chain(re_chill)
            .chain(finalize);
        // debug!("PARAMS\tEPOCH\tRESIDUAL\tCONSISTENCY\tINTEGER");
        for (epoch, (alpha, beta, gamma)) in schedule.enumerate() {
            // debug!("PARAMS\t{}\t{}\t{}\t{}", epoch, alpha, beta, gamma);
            self.gradient.iter_mut().for_each(|x| *x = 0f64);
            self.momentum.iter_mut().for_each(|x| *x = 0f64);
            // TODO: this is the worst hack I've ever implemented. Shame on me.
            copy_num
                .iter_mut()
                .for_each(|x| *x += (epoch as f64 + 1f64).recip());
            let mut loss_log: VecDeque<_> =
                std::iter::repeat(std::f64::INFINITY).take(100).collect();
            loop {
                let total_loss = self.update(copy_num, alpha, beta, gamma);
                let old = loss_log.pop_front().unwrap();
                loss_log.push_back(total_loss);
                if old - total_loss < 0.0001 {
                    break;
                };
            }
        }
    }
    // Momentum method
    fn update(&mut self, copy_num: &mut [f64], alpha: f64, beta: f64, gamma: f64) -> f64 {
        let total_loss = self.update_gradient(copy_num, alpha, beta, gamma);
        let (learn_rate, moment_coef) = (0.01, 0f64);
        for ((x, d), moment) in copy_num
            .iter_mut()
            .zip(self.gradient.iter())
            .zip(self.momentum.iter_mut())
        {
            let diff = moment_coef * *moment - learn_rate * d;
            *x += diff;
            *moment = diff;
        }
        total_loss
    }
    fn update_gradient(&mut self, copy_num: &[f64], alpha: f64, beta: f64, gamma: f64) -> f64 {
        let mut node_residual = vec![0f64; self.nodes.len()];
        let mut node_penalty = vec![0f64; self.nodes.len()];
        for (&(from, fplus, to, tplus, _), cp) in self.edges.iter().zip(copy_num.iter()) {
            let coef = if fplus { -1f64 } else { 1f64 };
            node_penalty[from] += cp * coef;
            let coef = if tplus { -1f64 } else { 1f64 };
            node_penalty[to] += cp * coef;
            node_residual[from] += cp / 2f64;
            node_residual[to] += cp / 2f64;
        }
        for (idx, _) in self.is_tip.iter().enumerate().filter(|(_, &x)| x) {
            node_residual[idx] *= 2f64;
            node_penalty[idx] = 0f64;
        }
        node_residual
            .iter_mut()
            .zip(self.nodes.iter())
            .for_each(|(x, y)| *x -= y);
        let integer_penalty: Vec<_> = copy_num
            .iter()
            .map(|x| x - x.round())
            .map(|x| if x.abs() < 0.5 - 0.00001 { x } else { 0f64 })
            .collect();
        for ((g, edge), cp) in self
            .gradient
            .iter_mut()
            .zip(self.edges.iter())
            .zip(copy_num.iter())
        {
            *g = cp - edge.4;
        }
        let edge_res: f64 = self.gradient.iter().map(|x| x * x).sum();
        let node_res: f64 = node_residual.iter().map(|x| x * x).sum();
        let node_pen: f64 = node_penalty.iter().map(|x| x * x).sum();
        let int_pen: f64 = integer_penalty.iter().map(|x| x * x).sum();
        self.gradient.iter_mut().for_each(|g| *g *= 2f64);
        for ((grad, &(from, fplus, to, tplus, _)), int_pen) in self
            .gradient
            .iter_mut()
            .zip(self.edges.iter())
            .zip(integer_penalty.iter())
        {
            // Grad is currently the residual of the edge.
            let residual = node_residual[from] + node_residual[to] + *grad;
            let from_coef = if fplus { -1f64 } else { 1f64 };
            let to_coef = if tplus { -1f64 } else { 1f64 };
            let node_consist = 2f64 * (from_coef * node_penalty[from] + to_coef * node_penalty[to]);
            *grad = alpha * residual + beta * node_consist + gamma * int_pen;
        }
        let total_loss = alpha * (edge_res + node_res) + beta * node_pen + gamma * int_pen;
        // debug!(
        //     "LOSS\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
        //     edge_res, node_res, node_pen, int_pen, total_loss
        // );
        total_loss
    }
}

use std::collections::HashMap;
/// Estimate copy number of GFA file.
/// Each segment record should have `cv:i:` samtag and each edge segment record
/// should have `cv:i:` tag and `ln:i:` tag.
/// The `cov` parameter is the haplotype coverage,
/// the `len` parameter is the average length of the raw reads,
/// and `unit_len` parameter is the length of the unit.
/// If the assembly graph is gapless, `len` and `unit_len` would be 0.
/// After estimation, the estimated copy number would be added to gfa as `cp:i:` tag.
pub fn estimate_copy_number_on_gfa(gfa: &mut gfa::GFA, cov: f64, lens: &[usize], unit_len: usize) {
    let node_index: HashMap<_, _> = gfa
        .iter()
        .filter_map(|record| match &record.content {
            gfa::Content::Seg(node) => Some(&node.sid),
            _ => None,
        })
        .enumerate()
        .map(|(idx, id)| (id.clone(), idx))
        .collect();
    let calibrator = CoverageCalibrator::new(lens);
    let old_cov = cov;
    let cov = calibrator.calib_f64(cov, unit_len);
    debug!("Coverage\t{}\t{}", old_cov, cov);
    let nodes: Vec<_> = gfa
        .iter()
        .filter_map(|record| {
            if let gfa::Content::Seg(_) = &record.content {
                let coverage: usize = record
                    .tags
                    .iter()
                    .find(|tag| tag.inner.starts_with("cv"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|x| x.parse().ok())?;
                let weight = calibrator.calib(coverage, unit_len) / cov;
                // debug!("DUMP\t{}\t{:.2}\tNode", coverage, weight);
                Some(weight)
            } else {
                None
            }
        })
        .collect();
    assert_eq!(nodes.len(), node_index.len());
    let edges: Vec<_> = gfa
        .iter()
        .filter_map(|record| match &record.content {
            gfa::Content::Edge(edge) => {
                let from = node_index[&edge.sid1.id];
                let from_plus = edge.beg1.pos == 0;
                let to = node_index[&edge.sid2.id];
                let to_plus = edge.beg2.pos == 0;
                let tags = &record.tags;
                let coverage: usize = tags
                    .iter()
                    .find(|x| x.inner.starts_with("cv"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|cov| cov.parse().ok())
                    .expect(&format!("{:?}", record.tags));
                let gap_len: isize = tags
                    .iter()
                    .find(|x| x.inner.starts_with("ln"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|cov| cov.parse().ok())
                    .expect(&format!("{:?}", record.tags));
                let gap_len = (gap_len + 2 * unit_len as isize).max(0) as usize;
                let weight = calibrator.calib(coverage, gap_len) / cov;
                // debug!("DUMP\t{}\t{:.2}\t{:.2}\tEdge", coverage, weight, gap_len);
                Some((from, from_plus, to, to_plus, weight))
            }
            _ => None,
        })
        .collect();
    let (node_cp, edge_cp) = estimate_copy_number(&nodes, &edges);
    assert_eq!(nodes.len(), node_cp.len());
    let nodes = gfa.iter_mut().filter(|record| match &record.content {
        gfa::Content::Seg(_) => true,
        _ => false,
    });
    for (record, cp) in nodes.zip(node_cp) {
        record.tags.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
    }
    assert_eq!(edges.len(), edge_cp.len());
    let edges = gfa.iter_mut().filter(|record| match &record.content {
        gfa::Content::Edge(_) => true,
        _ => false,
    });
    for (record, cp) in edges.zip(edge_cp) {
        record.tags.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
    }
}

/// Estimate copy number.
/// Return copy number of nodes, and copy number of edges.
pub fn estimate_copy_number(nodes: &[f64], edges: &[Edge]) -> (Vec<usize>, Vec<usize>) {
    let (is_tip, is_isolated): (Vec<_>, Vec<_>) = {
        // 2 * node_index + is_head.
        let mut degree = vec![0; nodes.len() * 2];
        for &(from, f_plus, to, t_plus, _) in edges.iter() {
            degree[2 * from + f_plus as usize] += 1;
            degree[2 * to + t_plus as usize] += 1;
        }
        // If either position is 0, it is a tip.
        degree
            .chunks_exact(2)
            .map(|w| (w[0] == 0 || w[1] == 0, w[0] == 0 && w[1] == 0))
            .unzip()
    };
    {
        let tip = is_tip.iter().filter(|&&x| x).count();
        let isolated = is_isolated.iter().filter(|&&x| x).count();
        let all = nodes.len();
        debug!("ISOLATED\tTIP\tNORMAL");
        debug!("{}\t{}\t{}", isolated, tip - isolated, all - tip);
    }
    let mut copy_num: Vec<_> = edges.iter().map(|f| f.4).collect();
    let mut optimizer = Optimizer::new(&nodes, &edges, &is_tip);
    optimizer.optimize(&mut copy_num);
    let mut node_cp = vec![0f64; nodes.len()];
    for (&(from, _, to, _, _), cp) in edges.iter().zip(copy_num.iter()) {
        node_cp[from] += cp / 2f64;
        node_cp[to] += cp / 2f64;
    }
    for (i, _) in is_tip.iter().enumerate().filter(|(_, &x)| x) {
        node_cp[i] *= 2f64;
    }
    for (i, _) in is_isolated.iter().enumerate().filter(|(_, &x)| x) {
        node_cp[i] = nodes[i];
    }
    // for (i, cp) in copy_num.iter().enumerate() {
    //     debug!("ROUND\t{}\t{:.2}\t{}\tEDGE", edges[i].4, cp, cp.round());
    // }
    // for (i, cp) in node_cp.iter().enumerate() {
    //     debug!("ROUND\t{}\t{:.2}\t{}\tNODE", nodes[i], cp, cp.round());
    // }
    let copy_num: Vec<_> = copy_num.iter().map(|&x| x.round() as usize).collect();
    let node_cp: Vec<_> = node_cp.iter().map(|&x| x.round() as usize).collect();
    (node_cp, copy_num)
}

// /// Optimization by using OpEn.
// pub fn estimate_copy_number_open(nodes: &[f64], edges: &[Edge]) -> (Vec<usize>, Vec<usize>) {
//     let (is_tip, is_isolated): (Vec<_>, Vec<_>) = {
//         // 2 * node_index + is_head.
//         let mut degree = vec![0; nodes.len() * 2];
//         for &(from, f_plus, to, t_plus, _) in edges.iter() {
//             degree[2 * from + f_plus as usize] += 1;
//             degree[2 * to + t_plus as usize] += 1;
//         }
//         // If either position is 0, it is a tip.
//         degree
//             .chunks_exact(2)
//             .map(|w| (w[0] == 0 || w[1] == 0, w[0] == 0 && w[1] == 0))
//             .unzip()
//     };
//     let (alpha, beta, gamma) = (0.1, 0.5, 0.5);
//     let loss_function =
//         |copy_num: &[f64], cost: &mut f64| -> Result<(), optimization_engine::SolverError> {
//             let mut node_residual: Vec<_> = nodes.iter().map(|x| -x).collect();
//             let mut node_penalty = vec![0f64; nodes.len()];
//             for (&(from, fplus, to, tplus, _), cp) in edges.iter().zip(copy_num.iter()) {
//                 let coef = if fplus { -1f64 } else { 1f64 };
//                 node_penalty[from] += cp * coef;
//                 let coef = if tplus { -1f64 } else { 1f64 };
//                 node_penalty[to] += cp * coef;
//                 node_residual[from] += cp / 2f64;
//                 node_residual[to] += cp / 2f64;
//             }
//             for (idx, _) in is_tip.iter().enumerate().filter(|(_, &x)| x) {
//                 node_residual[idx] *= 2f64;
//                 node_penalty[idx] = 0f64;
//             }
//             let node_res: f64 = node_residual.iter().map(|x| x * x).sum();
//             let node_pen: f64 = node_penalty.iter().map(|x| x * x).sum();

//             let integer_penalty = copy_num.iter().map(|x| x - x.round());
//             let int_pen: f64 = integer_penalty.map(|x| x * x).sum();
//             let edge_residual = edges.iter().zip(copy_num.iter()).map(|(e, c)| e.4 - c);
//             let edge_res: f64 = edge_residual.map(|x| x * x).sum();
//             *cost = alpha * (edge_res + node_res) + beta * node_pen + gamma * int_pen;
//             Ok(())
//         };
//     let gradient = |copy_num: &[f64],
//                     grad: &mut [f64]|
//      -> Result<(), optimization_engine::SolverError> {
//         let mut node_residual = vec![0f64; nodes.len()];
//         let mut node_penalty = vec![0f64; nodes.len()];
//         for (&(from, fplus, to, tplus, _), cp) in edges.iter().zip(copy_num.iter()) {
//             let coef = if fplus { -1f64 } else { 1f64 };
//             node_penalty[from] += cp * coef;
//             let coef = if tplus { -1f64 } else { 1f64 };
//             node_penalty[to] += cp * coef;
//             node_residual[from] += cp / 2f64;
//             node_residual[to] += cp / 2f64;
//         }
//         for (idx, _) in is_tip.iter().enumerate().filter(|(_, &x)| x) {
//             node_residual[idx] *= 2f64;
//             node_penalty[idx] = 0f64;
//         }
//         node_residual
//             .iter_mut()
//             .zip(nodes.iter())
//             .for_each(|(x, y)| *x -= y);
//         let integer_penalty: Vec<_> = copy_num
//             .iter()
//             .map(|x| x - x.round())
//             .map(|x| if x.abs() < 0.5 - 0.00001 { x } else { 0f64 })
//             .collect();
//         for ((g, edge), cp) in grad.iter_mut().zip(edges.iter()).zip(copy_num.iter()) {
//             *g = cp - edge.4;
//         }
//         grad.iter_mut().for_each(|g| *g *= 2f64);
//         for ((grad, &(from, fplus, to, tplus, _)), int_pen) in grad
//             .iter_mut()
//             .zip(edges.iter())
//             .zip(integer_penalty.iter())
//         {
//             // Grad is currently the residual of the edge.
//             let residual = node_residual[from] + node_residual[to] + *grad;
//             let from_coef = if fplus { -1f64 } else { 1f64 };
//             let to_coef = if tplus { -1f64 } else { 1f64 };
//             let node_consist = 2f64 * (from_coef * node_penalty[from] + to_coef * node_penalty[to]);
//             *grad = alpha * residual + beta * node_consist + gamma * int_pen;
//         }
//         Ok(())
//     };
//     use optimization_engine::Optimizer;
//     let mut copy_num: Vec<_> = edges.iter().map(|f| f.4).collect();
//     let constraint = optimization_engine::constraints::NoConstraints::new();
//     let problem = optimization_engine::Problem::new(&constraint, gradient, loss_function);
//     let mut cache = optimization_engine::core::panoc::PANOCCache::new(copy_num.len(), 1e-1, 10);
//     let mut optim = optimization_engine::core::panoc::PANOCOptimizer::new(problem, &mut cache);
//     debug!("{:?}", optim.solve(&mut copy_num).unwrap());
//     let mut loss = 0f64;
//     loss_function(&copy_num, &mut loss).unwrap();
//     debug!("Final loss\t{}", loss);
//     let mut node_cp = vec![0f64; nodes.len()];
//     for (&(from, _, to, _, _), cp) in edges.iter().zip(copy_num.iter()) {
//         node_cp[from] += cp / 2f64;
//         node_cp[to] += cp / 2f64;
//     }
//     for (i, _) in is_tip.iter().enumerate().filter(|(_, &x)| x) {
//         node_cp[i] *= 2f64;
//     }
//     for (i, _) in is_isolated.iter().enumerate().filter(|(_, &x)| x) {
//         node_cp[i] = nodes[i];
//     }
//     let copy_num: Vec<_> = copy_num.iter().map(|&x| x.round() as usize).collect();
//     let node_cp: Vec<_> = node_cp.iter().map(|&x| x.round() as usize).collect();
//     (node_cp, copy_num)
// }

#[cfg(test)]
mod cpe_test {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
    fn dataset1() -> (Vec<f64>, Vec<Edge>) {
        let nodes = vec![2f64, 1f64, 1f64, 2f64];
        let edges = vec![
            (0, false, 1, true, 1f64),
            (0, false, 2, true, 1f64),
            (1, false, 3, true, 1f64),
            (2, false, 3, true, 1f64),
        ];
        (nodes, edges)
    }
    #[test]
    fn dataset1_test() {
        let (nodes, edges) = dataset1();
        let (nodes_cp, edges_cp) = estimate_copy_number(&nodes, &edges);
        assert_eq!(nodes_cp, vec![2, 1, 1, 2]);
        assert_eq!(edges_cp, vec![1, 1, 1, 1]);
    }
    fn dataset2() -> (Vec<f64>, Vec<usize>, Vec<Edge>, Vec<usize>) {
        let nodes = vec![1f64, 1f64, 1f64, 1f64];
        let nodes_answer = vec![1; 4];
        let edges = vec![
            (0, false, 1, true, 1f64),
            (1, false, 2, true, 1f64),
            (2, false, 3, true, 1f64),
            (3, false, 0, true, 1f64),
        ];
        let edges_answer = vec![1; 4];
        (nodes, nodes_answer, edges, edges_answer)
    }
    #[test]
    fn dataset2_test() {
        let (n, na, e, ea) = dataset2();
        let (nc, ec) = estimate_copy_number(&n, &e);
        assert_eq!(nc, na);
        assert_eq!(ec, ea);
    }
    fn dataset3() -> (Vec<f64>, Vec<usize>, Vec<Edge>, Vec<usize>) {
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(4324);
        let nodes_answer = vec![2, 1, 2, 1, 1, 2, 4, 2];
        let nodes: Vec<_> = nodes_answer
            .iter()
            .map(|&cp| {
                let diff = rng.gen_range(-0.3..0.3);
                cp as f64 + diff
            })
            .collect();
        let edge_answer = vec![1, 1, 1, 1, 1, 1, 1, 2, 2, 2];
        let edges: Vec<_> = vec![
            (0, false, 1, true, 1),
            (0, false, 2, true, 1),
            (1, false, 2, true, 1),
            (2, false, 3, true, 1),
            (2, false, 4, true, 1),
            (3, false, 5, true, 1),
            (4, false, 5, true, 1),
            (5, false, 6, true, 2),
            (6, false, 6, true, 2),
            (6, false, 7, true, 2),
        ]
        .into_iter()
        .map(|(from, f, to, t, cov)| {
            let diff = rng.gen_range(-0.3..0.3);
            (from, f, to, t, cov as f64 + diff)
        })
        .collect();
        (nodes, nodes_answer, edges, edge_answer)
    }
    #[test]
    fn dataset3_test() {
        let (n, na, e, ea) = dataset3();
        let (nc, ec) = estimate_copy_number(&n, &e);
        assert_eq!(nc, na);
        assert_eq!(ec, ea);
    }
}

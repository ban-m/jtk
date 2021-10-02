use crate::assemble::{ditch_graph::DitchGraph, *};
use definitions::DataSet;
use serde::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiplicityEstimationConfig {
    max_cluster: usize,
    seed: u64,
    path: Option<String>,
    thread: usize,
}

impl MultiplicityEstimationConfig {
    pub fn new(thread: usize, max_cluster: usize, seed: u64, path: Option<&str>) -> Self {
        Self {
            thread,
            max_cluster,
            seed,
            path: path.map(|x| x.to_string()),
        }
    }
}

pub trait MultiplicityEstimation {
    // fn estimate_multiplicity(self, config: &MultiplicityEstimationConfig) -> Self;
    fn estimate_multiplicity_graph(self, config: &MultiplicityEstimationConfig) -> Self;
}

impl MultiplicityEstimation for DataSet {
    // fn estimate_multiplicity(mut self, config: &MultiplicityEstimationConfig) -> Self {
    //     for read in self.encoded_reads.iter_mut() {
    //         for node in read.nodes.iter_mut() {
    //             node.cluster = 0;
    //         }
    //     }
    //     self.assignments = self
    //         .encoded_reads
    //         .iter()
    //         .map(|r| definitions::Assignment::new(r.id, 0))
    //         .collect();
    //     let assemble_config = super::AssembleConfig::new(config.thread, 100, false, false);
    //     debug!("Start assembling {} reads", self.encoded_reads.len());
    //     let graph = self.assemble_draft_graph(&assemble_config);
    //     debug!("Assembled reads.");
    //     if let Some(mut file) = config
    //         .path
    //         .as_ref()
    //         .and_then(|path| std::fs::File::create(path).ok())
    //         .map(std::io::BufWriter::new)
    //     {
    //         use std::io::Write;
    //         let gfa = self.assemble(&assemble_config);
    //         writeln!(&mut file, "{}", gfa).unwrap();
    //     }
    //     debug!("GRAPH\tID\tCoverage\tMean\tLen");
    //     let estimated_cluster_num: HashMap<u64, usize> = {
    //         let (result, single_copy_coverage) = estimate_graph_multiplicity(&self, &graph, config);
    //         self.coverage = Some(single_copy_coverage);
    //         let mut cluster_num = HashMap::new();
    //         for (unit, cluster) in result {
    //             cluster_num.insert(unit, cluster);
    //         }
    //         cluster_num
    //     };
    //     for unit in self.selected_chunks.iter_mut() {
    //         if let Some(&cl_num) = estimated_cluster_num.get(&unit.id) {
    //             unit.cluster_num = cl_num;
    //         }
    //     }
    //     self
    // }
    fn estimate_multiplicity_graph(mut self, config: &MultiplicityEstimationConfig) -> Self {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *counts.entry(node.unit).or_default() += 1;
        }
        let cov = {
            let mut counts: Vec<_> = counts.values().copied().collect();
            counts.sort_unstable();
            counts[counts.len() / 2] as f64 / 2f64
        };
        debug!("MULTP\tCOVERAGE\t{}\tHAPLOID", cov);
        self.coverage = Some(cov);
        let reads: Vec<_> = self.encoded_reads.iter().collect();
        let assemble_config = AssembleConfig::new(config.thread, 100, false, false, 6);
        let mut graph =
            ditch_graph::DitchGraph::new(&reads, Some(&self.selected_chunks), &assemble_config);
        graph.remove_lightweight_edges(3, true);
        let lens: Vec<_> = self.raw_reads.iter().map(|x| x.seq().len()).collect();
        graph.assign_copy_number(cov, &lens);
        graph.remove_zero_copy_elements(&lens, 0.3);
        let (nodes, _) = graph.copy_number_estimation(cov, &lens);
        for chunk in self.selected_chunks.iter_mut() {
            chunk.cluster_num = match nodes.get(&(chunk.id, 0)) {
                Some(res) => *res,
                None => {
                    warn!("MULTP\tCHUNK_ID\t{}\tMissing\tFilled with 1", chunk.id);
                    1
                }
            };
        }
        let mut counts_group: HashMap<_, Vec<_>> = HashMap::new();
        for ((node, _), cp) in nodes.iter() {
            let occ = counts.get(node).unwrap_or(&0);
            counts_group.entry(cp).or_default().push(*occ as usize);
        }
        let mut counts_group: Vec<_> = counts_group.into_iter().collect();
        counts_group.sort_by_key(|x| x.0);
        for (cp, occs) in counts_group {
            let sum: usize = occs.iter().sum();
            let sumsq: usize = occs.iter().map(|x| x * x).sum();
            let mean = sum as f64 / occs.len() as f64;
            let sd = (sumsq as f64 / occs.len() as f64 - mean * mean).sqrt();
            debug!("MULTP\t{}\t{}\t{}\t{}", cp, occs.len(), mean, sd);
        }
        if let Some(mut file) = config
            .path
            .as_ref()
            .and_then(|path| std::fs::File::create(path).ok())
            .map(std::io::BufWriter::new)
        {
            use std::io::Write;
            let gfa = convert_to_gfa(&graph, &assemble_config);
            writeln!(&mut file, "{}", gfa).unwrap();
        }
        self
        // self.assignments = self
        //     .encoded_reads
        //     .iter()
        //     .map(|r| definitions::Assignment::new(r.id, 0))
        //     .collect();
        // self.encoded_reads
        //     .iter_mut()
        //     .for_each(|read| read.nodes.iter_mut().for_each(|n| n.cluster = 0));
        // debug!("Start assembling {} reads", self.encoded_reads.len());
        // let (mut gfa, tigname) = assemble_with_tigname(&self, config);
        // let mut counts: HashMap<_, u32> = HashMap::new();
        // for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        //     *counts.entry(node.unit).or_default() += 1;
        // }
        // let cov = {
        //     let mut counts: Vec<_> = counts.values().copied().collect();
        //     counts.sort_unstable();
        //     counts[counts.len() / 2] as f64 / 2f64
        // };
        // self.coverage = Some(cov);
        // let (unit_len, unit_num) = self
        //     .selected_chunks
        //     .iter()
        //     .fold((0, 0), |(tot, num), u| (tot + u.seq().len(), num + 1));
        // let unit_len = unit_len / unit_num;
        // let lens: Vec<_> = self
        //     .encoded_reads
        //     .iter()
        //     .map(|r| r.original_length)
        //     .collect();
        // debug!("HAPLOID\t{}", cov);
        // debug!("ESTIM\tTIGID\tLEN\tCP\tCOV");
        // crate::assemble::copy_number::estimate_copy_number_on_gfa(&mut gfa, cov, &lens, unit_len);
        // let mut estimated_cluster_num: HashMap<_, _> = HashMap::new();
        // for record in gfa.iter() {
        //     if let gfa::Content::Seg(seg) = &record.content {
        //         let cp: usize = record
        //             .tags
        //             .iter()
        //             .find(|x| x.inner.starts_with("cp"))
        //             .and_then(|tag| tag.inner.split(':').nth(2))
        //             .and_then(|cp| cp.parse().ok())
        //             .unwrap_or_else(|| panic!("{:?}", record.tags));
        //         let units = tigname.get(&seg.sid).unwrap();
        //         let cov: u32 = units.iter().map(|u| counts[u]).sum();
        //         let cov = cov / units.len() as u32;
        //         debug!("ESTIM\t{}\t{}\t{}\t{}", seg.sid, units.len(), cp, cov);
        //         for &unit in units.iter() {
        //             estimated_cluster_num.insert(unit, cp.max(1));
        //         }
        //     }
        // }
        // for unit in self.selected_chunks.iter_mut() {
        //     if let Some(&cl_num) = estimated_cluster_num.get(&unit.id) {
        //         unit.cluster_num = cl_num;
        //     }
        // }
        // if let Some(mut file) = config
        //     .path
        //     .as_ref()
        //     .and_then(|path| std::fs::File::create(path).ok())
        //     .map(std::io::BufWriter::new)
        // {
        //     use std::io::Write;
        //     writeln!(&mut file, "{}", gfa).unwrap();
        // }
        // self
    }
}

fn convert_to_gfa(graph: &DitchGraph, c: &AssembleConfig) -> gfa::GFA {
    let (segments, edge, group, summaries) = graph.spell(c, 0);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("{} segments({} bp in total).", segments.len(), total_base);
    // TODO: maybe just zip up segments and summaries would be OK?
    let nodes = segments.into_iter().map(|node| {
        let tags = summaries
            .iter()
            .find(|x| x.id == node.sid)
            .map(|contigsummary| {
                let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
                let coverage =
                    gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
                let (cp, cpnum) = contigsummary
                    .summary
                    .iter()
                    .filter_map(|elm| elm.copy_number)
                    .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
                let mut tags = vec![coverage];
                if cpnum != 0 {
                    tags.push(gfa::SamTag::new(format!("cp:i:{}", cp / cpnum)));
                }
                tags
            })
            .unwrap_or_else(Vec::new);
        gfa::Record::from_contents(gfa::Content::Seg(node), tags)
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    let group = std::iter::once(group);
    let header = gfa::Content::Header(gfa::Header::default());
    let header = std::iter::once(gfa::Record::from_contents(header, vec![]));
    let records: Vec<_> = header.chain(group).chain(nodes).chain(edges).collect();
    gfa::GFA::from_records(records)
}

// fn assemble_with_tigname(
//     ds: &DataSet,
//     config: &MultiplicityEstimationConfig,
// ) -> (gfa::GFA, HashMap<String, Vec<u64>>) {
//     let assemble_config = super::AssembleConfig::new(config.thread, 100, false, false, 6);
//     debug!("Start assembly");
//     let header = gfa::Content::Header(gfa::Header::default());
//     let header = gfa::Record::from_contents(header, vec![]);
//     use crate::assemble::assemble;
//     let (records, summaries) = assemble(ds, 0, &assemble_config);
//     let mut header = vec![header];
//     header.extend(records);
//     let tigname: HashMap<_, _> = summaries
//         .into_iter()
//         .map(
//             |crate::assemble::ditch_graph::ContigSummary { id, summary }| {
//                 let summary: Vec<_> = summary.iter().map(|s| s.unit).collect();
//                 (id, summary)
//             },
//         )
//         .collect();
//     (gfa::GFA::from_records(header), tigname)
// }

// fn estimate_graph_multiplicity(
//     ds: &DataSet,
//     graph: &super::assemble::Graph,
//     c: &MultiplicityEstimationConfig,
// ) -> (Vec<(u64, usize)>, f64) {
//     let covs: Vec<_> = graph
//         .nodes
//         .iter()
//         .map(|node| {
//             let len = node.segments.len();
//             let unit: HashSet<_> = node.segments.iter().map(|t| t.unit).collect();
//             let coverage = ds
//                 .encoded_reads
//                 .iter()
//                 .map(|r| r.nodes.iter().filter(|n| unit.contains(&n.unit)).count())
//                 .sum::<usize>();
//             let mean = (coverage / len) as u64;
//             mean
//         })
//         .collect();
//     use rayon::prelude::*;
//     let (model, aic): (Model, f64) = (1..c.max_cluster)
//         .into_par_iter()
//         .map(|k| {
//             let seed = k as u64 + c.seed;
//             let (model, lk) = clustering(&covs, k, seed);
//             // Lambda for each cluster, fraction for each cluster,
//             // and one constraint that the sum of the fractions equals to 1.
//             let aic = -2. * lk + (2. * k as f64 - 1.);
//             (model, aic)
//         })
//         .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
//         .unwrap();
//     debug!("AIC\t{}", aic);
//     let assignments: Vec<_> = covs.iter().map(|&d| model.assign(d)).collect();
//     let single_copy_coverage = {
//         let min = model
//             .lambdas
//             .iter()
//             .min_by(|a, b| a.partial_cmp(&b).unwrap())
//             .unwrap();
//         let coverage = model
//             .lambdas
//             .iter()
//             .filter(|&x| (x - min).abs() > 0.01)
//             .min_by(|a, b| a.partial_cmp(&b).unwrap())
//             .unwrap_or(min);
//         debug!("DIPLOID\t{}", coverage);
//         coverage / 2.
//     };
//     let repeat_num: Vec<_> = model
//         .lambdas
//         .iter()
//         .map(|x| ((x / single_copy_coverage) + 0.5).floor() as usize)
//         .collect();
//     debug!("LAMBDAS:{:?}", model.lambdas);
//     debug!("PREDCT:{:?}", repeat_num);
//     let mut result = vec![];
//     debug!("REPEATNUM\tID\tMULTP\tCLUSTER");
//     for (&cl, contig) in assignments.iter().zip(graph.nodes.iter()) {
//         let repeat_num = repeat_num[cl];
//         debug!("REPEATNUM\t{}\t{}\t{}", contig.id, repeat_num, cl);
//         for node in contig.segments.iter() {
//             result.push((node.unit, repeat_num));
//         }
//     }
//     (result, single_copy_coverage)
// }

// struct Model {
//     cluster: usize,
//     fractions: Vec<f64>,
//     lambdas: Vec<f64>,
// }
// const SMALL: f64 = 0.00000000000000001;
// impl Model {
//     fn new(data: &[u64], weight: &[Vec<f64>], k: usize) -> Self {
//         let sum: Vec<_> = (0..k)
//             .map(|cl| weight.iter().map(|ws| ws[cl]).sum::<f64>() + SMALL)
//             .collect();
//         let fractions: Vec<_> = sum.iter().map(|w| w / data.len() as f64).collect();
//         let lambdas: Vec<_> = sum
//             .iter()
//             .enumerate()
//             .map(|(cl, sum)| {
//                 weight
//                     .iter()
//                     .zip(data)
//                     .map(|(ws, &x)| x as f64 * ws[cl])
//                     .sum::<f64>()
//                     / sum
//             })
//             .collect();
//         let cluster = k;
//         Self {
//             cluster,
//             fractions,
//             lambdas,
//         }
//     }
//     fn lk(&self, data: &[u64]) -> f64 {
//         data.iter().map(|&d| self.lk_data(d)).sum::<f64>()
//     }
//     fn lk_data(&self, data: u64) -> f64 {
//         let lks: Vec<_> = (0..self.cluster)
//             .map(|cl| {
//                 self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//                     - self.lambdas[cl]
//                     - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
//             })
//             .collect();
//         logsumexp(&lks)
//     }
//     fn update_weight(&self, ws: &mut [f64], data: u64) {
//         assert_eq!(ws.len(), self.cluster);
//         for (cl, w) in ws.iter_mut().enumerate() {
//             *w = self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//                 - self.lambdas[cl]
//                 - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>();
//         }
//         let lk = logsumexp(&ws);
//         ws.iter_mut().for_each(|x| *x = (*x - lk).exp());
//     }
//     // fn new_weight(&self, data: u64) -> Vec<f64> {
//     //     let lks: Vec<_> = (0..self.cluster)
//     //         .map(|cl| {
//     //             self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//     //                 - self.lambdas[cl]
//     //                 - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
//     //         })
//     //         .collect();
//     //     let lk = logsumexp(&lks);
//     //     assert!((1. - lks.iter().map(|x| (x - lk).exp()).sum::<f64>()).abs() < 0.0001);
//     //     lks.iter().map(|x| (x - lk).exp()).collect()
//     // }
//     fn assign(&self, data: u64) -> usize {
//         let (cl, _) = (0..self.cluster)
//             .map(|cl| {
//                 self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//                     - self.lambdas[cl]
//                     - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
//             })
//             .enumerate()
//             .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
//             .unwrap();
//         cl
//     }
// }

// fn logsumexp(xs: &[f64]) -> f64 {
//     if xs.is_empty() {
//         return 0.;
//     }
//     let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
//     let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
//     assert!(sum >= 0., "{:?}->{}", xs, sum);
//     max + sum
// }

// fn clustering(data: &[u64], k: usize, seed: u64) -> (Model, f64) {
//     let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
//     let (weight, lk) = (0..5)
//         .map(|_| {
//             let mut weight: Vec<_> = (0..data.len())
//                 .map(|_| {
//                     let mut ws = vec![0.; k];
//                     ws[rng.gen::<usize>() % k] = 1.;
//                     ws
//                 })
//                 .collect();
//             let mut lk = std::f64::NEG_INFINITY;
//             loop {
//                 let model = Model::new(data, &weight, k);
//                 let new_lk = model.lk(data);
//                 let diff = new_lk - lk;
//                 if diff < 0.00001 {
//                     break;
//                 }
//                 lk = new_lk;
//                 for (ws, d) in weight.iter_mut().zip(data.iter()) {
//                     model.update_weight(ws, *d);
//                 }
//             }
//             (weight, lk)
//         })
//         .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
//         .unwrap();
//     let model = Model::new(data, &weight, k);
//     (model, lk)
// }

// // pub fn cluster_coverage(
// //     unit_covs: &HashMap<u64, u64>,
// //     c: &MultiplicityEstimationConfig,
// // ) -> (Vec<(u64, usize)>, f64) {
// //     use rayon::prelude::*;
// //     let covs: Vec<_> = unit_covs.values().copied().collect();
// //     let (model, aic): (Model, f64) = (1..c.max_cluster)
// //         .into_par_iter()
// //         .map(|k| {
// //             let seed = k as u64 + c.seed;
// //             let (model, lk) = clustering(&covs, k, seed);
// //             let aic = -2. * lk + (2. * k as f64 - 1.);
// //             (model, aic)
// //         })
// //         .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
// //         .unwrap();
// //     debug!("AIC\t{}", aic);
// //     let assignments: Vec<_> = covs.iter().map(|&d| model.assign(d)).collect();
// //     let single_copy_coverage = {
// //         let min = model
// //             .lambdas
// //             .iter()
// //             .min_by(|a, b| a.partial_cmp(&b).unwrap())
// //             .unwrap();
// //         let coverage = model
// //             .lambdas
// //             .iter()
// //             .filter(|&x| (x - min).abs() > 0.01)
// //             .min_by(|a, b| a.partial_cmp(&b).unwrap())
// //             .unwrap_or(min);
// //         debug!("DIPLOID\t{}", coverage);
// //         coverage / 2.
// //     };
// //     let repeat_num: Vec<_> = model
// //         .lambdas
// //         .iter()
// //         .map(|x| ((x / single_copy_coverage) + 0.5).floor() as usize)
// //         .collect();
// //     debug!("LAMBDAS:{:?}", model.lambdas);
// //     debug!("PREDCT:{:?}", repeat_num);
// //     let mut result = vec![];
// //     for (&cl, (&unit, &_)) in assignments.iter().zip(unit_covs.iter()) {
// //         let repeat_num = repeat_num[cl];
// //         result.push((unit, repeat_num));
// //     }
// //     (result, single_copy_coverage)
// // }

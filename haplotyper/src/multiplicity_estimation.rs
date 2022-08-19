use crate::assemble::{ditch_graph::DitchGraph, *};
use definitions::DataSet;
use serde::*;
use std::collections::HashMap;
use std::collections::HashSet;
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiplicityEstimationConfig {
    seed: u64,
    path: Option<String>,
}

impl MultiplicityEstimationConfig {
    pub fn new(seed: u64, path: Option<&str>) -> Self {
        Self {
            seed,
            path: path.map(|x| x.to_string()),
        }
    }
}

pub trait MultiplicityEstimation {
    fn estimate_multiplicity(&mut self, config: &MultiplicityEstimationConfig);
    // Remove units with copy number more than or equal to `upper`
    fn purge_multiplicity(&mut self, upper: usize);
}

impl MultiplicityEstimation for DataSet {
    fn estimate_multiplicity(&mut self, config: &MultiplicityEstimationConfig) {
        crate::misc::update_coverage(self);
        let cov = self.coverage.unwrap();
        let reads: Vec<_> = self.encoded_reads.iter().collect();
        let assemble_config = AssembleConfig::new(100, false, false, 4, 0f64, false);
        let rt = self.read_type;
        let mut graph =
            ditch_graph::DitchGraph::new(&reads, &self.selected_chunks, rt, &assemble_config);
        let thr = match self.read_type {
            definitions::ReadType::CCS => 1,
            definitions::ReadType::CLR => 2,
            definitions::ReadType::ONT => 2,
            definitions::ReadType::None => 1,
        };
        graph.remove_lightweight_edges(thr, true);
        debug!("SQUISHED\t{graph}");
        use rand::SeedableRng;
        use rand_xoshiro::Xoroshiro128PlusPlus;
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(config.seed);
        // graph.assign_copy_number_mcmc(cov, &mut rng);
        // graph.assign_copy_number_mst(cov, &mut rng);
        graph.assign_copy_number_flow(cov, &mut rng);
        let nodes: HashMap<_, _> = graph
            .nodes()
            .filter_map(|(_, node)| node.copy_number.map(|c| (node.node, c)))
            .collect();
        let mut chunks: HashMap<_, _> =
            self.selected_chunks.iter_mut().map(|c| (c.id, c)).collect();
        // reset copy number.
        chunks.values_mut().for_each(|c| c.copy_num = 0);
        for ((chunk, _), cp) in nodes.iter() {
            if let Some(c) = chunks.get_mut(chunk) {
                c.copy_num += cp;
            }
        }
        // Rescue 0-copy number chunks.
        chunks
            .values_mut()
            .for_each(|c| c.copy_num = c.copy_num.max(1));
        {
            // DUMP Information
            let mut counts: HashMap<(u64, u64), usize> = HashMap::new();
            for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
                *counts.entry((node.unit, node.cluster)).or_default() += 1;
            }
            let mut counts_group: HashMap<_, Vec<_>> = HashMap::new();
            for (node, cp) in nodes.iter() {
                let occ = counts.get(node).unwrap_or(&0);
                counts_group.entry(cp).or_default().push(*occ as usize);
            }
            let mut counts_group: Vec<_> = counts_group.into_iter().collect();
            counts_group.sort_by_key(|x| x.0);
            debug!("MULTP\tCopyNum\tOccs\tMean\tSD");
            for (cp, occs) in counts_group {
                let sum: usize = occs.iter().sum();
                let sumsq: usize = occs.iter().map(|x| x * x).sum();
                let mean = sum as f64 / occs.len() as f64;
                let sd = (sumsq as f64 / occs.len() as f64 - mean * mean).sqrt();
                debug!("MULTP\t{}\t{}\t{:.2}\t{:.2}", cp, occs.len(), mean, sd);
            }
        }
        if let Some(mut file) = config
            .path
            .as_ref()
            .and_then(|path| std::fs::File::create(path).ok())
            .map(std::io::BufWriter::new)
        {
            use std::io::Write;
            let gfa = convert_to_gfa(&graph, &assemble_config);
            writeln!(file, "{}", gfa).unwrap();
        }
    }
    fn purge_multiplicity(&mut self, upper: usize) {
        let to_remove: HashSet<_> = self
            .selected_chunks
            .iter()
            .filter_map(|c| (upper <= c.copy_num || c.copy_num == 0).then(|| c.id))
            .collect();
        {
            let mut counts: HashMap<_, u32> = HashMap::new();
            for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
                *counts.entry(node.unit).or_default() += 1;
            }
            for chunk in self
                .selected_chunks
                .iter()
                .filter(|c| to_remove.contains(&c.id))
            {
                let count = counts[&chunk.id];
                debug!("REMOVE\t{}\t{count}\t{}", chunk.id, chunk.copy_num);
            }
        }
        self.selected_chunks.retain(|x| !to_remove.contains(&x.id));
        use rayon::prelude::*;
        self.encoded_reads.par_iter_mut().for_each(|read| {
            let mut idx = 0;
            loop {
                match read.nodes.get(idx) {
                    Some(node) if to_remove.contains(&node.unit) => read.remove(idx),
                    Some(_) => idx += 1,
                    None => break,
                }
            }
        });
        self.encoded_reads.retain(|read| !read.nodes.is_empty());
    }
}

fn convert_to_gfa(graph: &DitchGraph, c: &AssembleConfig) -> gfa::GFA {
    let (segments, edge, _group, summaries, _) = graph.spell(c);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("MULTIP\tAssembly\t{}\t{}bp", segments.len(), total_base);
    let mut groups: HashMap<_, Vec<_>> = HashMap::new();
    let nodes: Vec<_> = segments
        .into_iter()
        .map(|node| {
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
                    groups
                        .entry(cp / cpnum.max(1))
                        .or_default()
                        .push(node.sid.clone());
                    tags
                })
                .unwrap_or_else(Vec::new);
            gfa::Record::from_contents(gfa::Content::Seg(node), tags.into())
        })
        .collect();
    {
        for summary in summaries.iter() {
            let (copy_num, tig_num) = summary
                .summary
                .iter()
                .filter_map(|s| s.copy_number)
                .fold((0, 0), |(c, x), copynum| (c + copynum, x + 1));
            let copy_num = match tig_num {
                0 => 0,
                _ => (copy_num as f64 / tig_num as f64).round() as usize,
            };

            let ids: Vec<_> = summary
                .summary
                .iter()
                .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
                .collect();
            debug!(
                "MULTIP\tContig\t{}\t{}\t{}",
                summary.id,
                copy_num,
                ids.join("\t")
            );
        }
    }
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags.into()));
    let groups = groups.into_iter().map(|(cp, ids)| {
        let group = gfa::UnorderedGroup {
            uid: Some(format!("cp:i:{}", cp)),
            ids,
        };
        let group = gfa::Content::Group(gfa::Group::Set(group));
        gfa::Record::from_contents(group, vec![].into())
    });
    // let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![].into());
    // let group = std::iter::once(group);
    let header = gfa::Content::Header(gfa::Header::default());
    let header = std::iter::once(gfa::Record::from_contents(header, vec![].into()));
    let records: Vec<_> = header.chain(groups).chain(nodes).chain(edges).collect();
    gfa::GFA::from_records(records)
}

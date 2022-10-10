pub mod copy_number;
pub mod ditch_graph;
use definitions::*;
use ditch_graph::*;
use gfa::GFA;
use serde::*;
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Graph {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node {
    pub id: String,
    pub segments: Vec<Tile>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tile {
    pub unit: u64,
    pub cluster: u64,
    pub strand: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    pub from: String,
    pub from_tail: bool,
    pub to: String,
    pub to_tail: bool,
}

impl Graph {
    pub fn enumerate_connected_components(&self) -> Vec<Vec<&Node>> {
        use crate::find_union::FindUnion;
        let mut fu = FindUnion::new(self.nodes.len());
        let id2index: HashMap<_, usize> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(index, node)| (&node.id, index))
            .collect();
        // Merge edges.
        for edge in self.edges.iter() {
            let from = id2index[&edge.from];
            let to = id2index[&edge.to];
            fu.unite(from, to);
        }
        // Take components.
        id2index
            .iter()
            .filter_map(|(_, &index)| {
                (fu.find(index).unwrap() == index).then(|| {
                    self.nodes
                        .iter()
                        .enumerate()
                        .filter_map(|(i, node)| (fu.find(i).unwrap() == index).then(|| node))
                        .collect::<Vec<&Node>>()
                })
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct AssembleConfig {
    to_polish: bool,
    window_size: usize,
    to_resolve: bool,
    min_span_reads: usize,
    span_likelihood_ratio: f64,
    to_bypass_contigs: bool,
}

impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {
            to_polish: false,
            window_size: 100,
            to_resolve: false,
            min_span_reads: 6,
            span_likelihood_ratio: 3f64,
            to_bypass_contigs: false,
        }
    }
}
impl AssembleConfig {
    pub fn new(
        window_size: usize,
        to_polish: bool,
        to_resolve: bool,
        min_span_reads: usize,
        span_likelihood_ratio: f64,
        to_bypass_contigs: bool,
    ) -> Self {
        Self {
            window_size,
            to_polish,
            to_resolve,
            min_span_reads,
            span_likelihood_ratio,
            to_bypass_contigs,
        }
    }
}

pub trait Assemble {
    /// Assemble the dataset. If there's duplicated regions or
    /// unresolved regions, it tries to un-entangle that region.
    fn assemble(&self, c: &AssembleConfig) -> GFA;
}

impl Assemble for DataSet {
    fn assemble(&self, c: &AssembleConfig) -> GFA {
        debug!("Start assembly");
        let (records, summaries) = assemble(self, c);
        let copy_numbers = get_contig_copy_numbers(&summaries);
        let shared_read_counts = count_contig_connection(self, &summaries);
        debug!("ContigConnection\tid1\tid2\tcp1\tcp2\tcount");
        for (i, cs) in shared_read_counts.iter().enumerate() {
            for (j, &c) in cs.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                let (iname, jname) = (&summaries[i].id, &summaries[j].id);
                let (icp, jcp) = (copy_numbers[i], copy_numbers[j]);
                debug!(
                    "ContigConnection\t{}\t{}\t{}\t{}\t{}",
                    iname, jname, icp, jcp, c
                );
            }
        }
        for summary in summaries {
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
            debug!("CONUNIT\t{}\t{}\t{}", summary.id, copy_num, ids.join("\t"));
        }
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![].into());
        let mut header = vec![header];
        header.extend(records);
        GFA::from_records(header)
    }
}

const LOWER_FRAC: f64 = 0.15;
/// ASSEMBLEIMPL
pub fn assemble(ds: &DataSet, c: &AssembleConfig) -> (Vec<gfa::Record>, Vec<ContigSummary>) {
    assert!(c.to_resolve);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let cov = ds.coverage.unwrap();
    let mut graph = DitchGraph::new(&reads, &ds.selected_chunks, ds.read_type, c);
    debug!("GRAPH\t{graph}");
    let thr = (cov * LOWER_FRAC).round() as usize;
    graph.remove_lightweight_edges(thr / 2 + 1, false);
    graph.remove_lightweight_edges(thr, true);
    // graph = {
    //     let mut old = graph.clone();
    //     old.remove_lightweight_edges(thr, true);
    //     graph.remove_lightweight_edges(thr, false);
    //     match graph.cc() {
    //         1 => graph,
    //         _ => old,
    //     }
    // };
    // graph.remove_lightweight_edges(thr, false);
    graph.clean_up_graph_for_assemble(cov, &reads, c, ds.read_type);
    let (mut segments, mut edges, _, summaries, encodings) = graph.spell(c);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("{} segments({} bp in total).", segments.len(), total_base);
    if c.to_polish {
        use crate::consensus;
        use crate::consensus::Polish;
        let seed = 394802;
        let radius = 50;
        let round = 4;
        let config =
            consensus::PolishConfig::new(seed, c.min_span_reads, c.window_size, radius, round);
        segments = ds.polish_segment(&segments, &encodings, &config);
        let lengths: HashMap<_, _> = segments
            .iter()
            .map(|seg| (seg.sid.clone(), seg.slen))
            .collect();
        for edge in edges.iter_mut() {
            if edge.0.beg1.pos != 0 {
                edge.0.beg1.pos = lengths[&edge.0.sid1.id] as usize;
                edge.0.end1.pos = lengths[&edge.0.sid1.id] as usize;
            }
            if edge.0.beg2.pos != 0 {
                edge.0.beg2.pos = lengths[&edge.0.sid2.id] as usize;
                edge.0.end2.pos = lengths[&edge.0.sid2.id] as usize;
            }
        }
    }
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
                    let copynum = (cp as f64 / cpnum.max(1) as f64).round() as usize;
                    let mut tags = vec![coverage];
                    if cpnum != 0 {
                        tags.push(gfa::SamTag::new(format!("cp:i:{copynum}")));
                    }
                    groups.entry(copynum).or_default().push(node.sid.clone());
                    tags
                })
                .unwrap_or_else(Vec::new);
            gfa::Record::from_contents(gfa::Content::Seg(node), tags.into())
        })
        .collect();
    let edges = edges
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
    let records: Vec<_> = groups.chain(nodes).chain(edges).collect();
    (records, summaries)
}

fn get_contig_copy_numbers(summaries: &[ContigSummary]) -> Vec<usize> {
    summaries
        .iter()
        .map(|summary| {
            let (cp, num) = summary
                .summary
                .iter()
                .filter_map(|n| n.copy_number)
                .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
            cp / num.max(1)
        })
        .collect()
}

/// Return the co-occurence of the reads.
/// [i][j] -> # of reads shared by summaries[i] and summaries[j].
fn count_contig_connection(ds: &DataSet, summaries: &[ContigSummary]) -> Vec<Vec<u32>> {
    // (unit,cluster) -> indices of the summary whose either terminal is (unit,cluster), if there is any.
    let mut contig_terminals: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for (i, summary) in summaries.iter().enumerate() {
        let mut summary = summary.summary.iter();
        let first = summary.next().unwrap();
        contig_terminals
            .entry((first.unit, first.cluster))
            .or_default()
            .push(i);
        if let Some(last) = summary.last() {
            contig_terminals
                .entry((last.unit, last.cluster))
                .or_default()
                .push(i);
        }
    }
    let mut shared_reads = vec![vec![0; summaries.len()]; summaries.len()];
    for read in ds.encoded_reads.iter() {
        let mut read_through = vec![];
        let mut is_first = true;
        for w in read.nodes.windows(2) {
            // If this (n,c) - (n',c') is connecting
            // two contigs, we record the both contigs.
            let (from, to) = match w {
                [f, t] => ((f.unit, f.cluster), (t.unit, t.cluster)),
                _ => unreachable!(),
            };
            if let (Some(f), Some(t)) = (contig_terminals.get(&from), contig_terminals.get(&to)) {
                // We record the leaving contig only at the first contig.
                // This is because, what we want to get is the list of arrived contig,
                // not the changelog of the contig.
                if is_first {
                    read_through.push(f);
                    is_first = false;
                }
                read_through.push(t);
            }
        }
        // If it is empty, it is NO-op, so it's ok.
        for &p1 in read_through.iter() {
            for &p2 in read_through.iter() {
                for &i in p1 {
                    for &j in p2 {
                        shared_reads[i][j] += 1;
                    }
                }
            }
        }
    }
    shared_reads
}

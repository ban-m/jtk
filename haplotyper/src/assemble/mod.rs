pub mod copy_number;
pub mod ditch_graph;
// pub mod string_graph;
// pub mod string_graph_dev;
use definitions::*;
use ditch_graph::*;
use gfa::GFA;
use serde::*;
use std::collections::{HashMap, HashSet};

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
    threads: usize,
    to_polish: bool,
    #[allow(dead_code)]
    window_size: usize,
    to_resolve: bool,
    min_span_reads: usize,
    span_likelihood_ratio: f64,
}

impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            to_polish: false,
            window_size: 100,
            to_resolve: false,
            min_span_reads: 6,
            span_likelihood_ratio: 3f64,
        }
    }
}
impl AssembleConfig {
    pub fn new(
        threads: usize,
        window_size: usize,
        to_polish: bool,
        to_resolve: bool,
        min_span_reads: usize,
        span_likelihood_ratio: f64,
    ) -> Self {
        Self {
            window_size,
            threads,
            to_polish,
            to_resolve,
            min_span_reads,
            span_likelihood_ratio,
        }
    }
}

pub trait Assemble {
    /// Assemble the dataset. If there's duplicated regions or
    /// unresolved regions, it tries to un-entangle that region.
    fn assemble(&self, c: &AssembleConfig) -> GFA;
    /// Assmeble the dataset into a draft assembly. It does not
    /// dive into difficult regions such as non-single-copy chunks or
    /// tangles. It just assemble the dataset into a graph.
    fn assemble_draft_graph(&self, c: &AssembleConfig) -> Graph;
    /// Detecting and squishing clusters with small confidence.
    fn squish_small_contig(&mut self, c: &AssembleConfig, len: usize);
    /// Zip up suspicious haplotig. In other words,
    /// after usual assembly workflow, if two haplotig shares reads more than `count`,
    /// it indicates that these contigs can not be phased because of something. We squish such haplotigs with length less than `len`.
    fn zip_up_suspicious_haplotig(&mut self, c: &AssembleConfig, count: u32, len: usize);
}

impl Assemble for DataSet {
    fn zip_up_suspicious_haplotig(&mut self, c: &AssembleConfig, count: u32, len: usize) {
        let (_, summaries) = assemble(self, c);
        // Select unique units.
        let copy_numbers = get_contig_copy_numbers(&summaries);
        let shared_read_counts = count_contig_connection(self, &summaries);
        // Squishing. Unit -> Vec<Cluster>
        let mut squishing_node: HashMap<u64, Vec<u64>> = HashMap::new();
        for (i, (cs, &i_copy_number)) in shared_read_counts
            .iter()
            .zip(copy_numbers.iter())
            .enumerate()
        {
            for (&c, &j_copy_number) in cs.iter().zip(copy_numbers.iter()) {
                if c < count || 2 <= j_copy_number || 2 <= i_copy_number {
                    continue;
                }
                if summaries[i].summary.len() < len {
                    let dump: Vec<_> = summaries[i]
                        .summary
                        .iter()
                        .map(|n| (n.unit, n.cluster))
                        .collect();
                    trace!("ZIPUP\t{:?}", dump);
                    for (u, c) in summaries[i].summary.iter().map(|n| (n.unit, n.cluster)) {
                        squishing_node.entry(u).or_default().push(c);
                    }
                    break;
                }
            }
        }
        // Squishing (unit,cluster)->cluster.
        let mut converting_map: HashMap<(u64, u64), u64> = HashMap::new();
        for (&unit, sqs) in squishing_node.iter() {
            let target: u64 = *sqs.iter().min().unwrap();
            if sqs.len() <= 2 {
                debug!("ZIPUP\t({},{:?})\t{}", unit, sqs, target);
                converting_map.extend(sqs.iter().map(|&cluster| ((unit, cluster), target)));
            }
        }
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if let Some(&to) = converting_map.get(&(node.unit, node.cluster)) {
                    node.cluster = to;
                }
            }
        }
        let reads: Vec<_> = self.encoded_reads.iter().collect();
        let lens: Vec<_> = self.raw_reads.iter().map(|x| x.seq().len()).collect();
        let cov = self.coverage.unwrap();
        let rt = self.read_type;
        let mut graph = DitchGraph::new(&reads, Some(&self.selected_chunks), rt, c);
        graph.remove_lightweight_edges(2, true);
        graph.clean_up_graph_for_assemble(cov, &lens, &reads, c, self.read_type);
        // TODO: Parametrize here.
        let squish = graph.squish_bubbles(2);
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .for_each(|n| {
                if let Some(res) = squish.get(&(n.unit, n.cluster)) {
                    n.cluster = *res;
                }
            });
    }
    fn squish_small_contig(&mut self, c: &AssembleConfig, len: usize) {
        // assert!(c.to_resolve);
        let reads: Vec<_> = self.encoded_reads.iter().collect();
        let cov = self.coverage.unwrap_or_else(|| panic!("Need coverage!"));
        let lens: Vec<_> = self.raw_reads.iter().map(|x| x.seq().len()).collect();
        let mut graph = DitchGraph::new(&reads, Some(&self.selected_chunks), self.read_type, c);
        match self.read_type {
            ReadType::CCS => graph.remove_lightweight_edges(1, true),
            ReadType::ONT | ReadType::None | ReadType::CLR => {
                graph.remove_lightweight_edges(2, true)
            }
        };
        graph.clean_up_graph_for_assemble(cov, &lens, &reads, c, self.read_type);
        let squish = graph.squish_bubbles(len);
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .for_each(|n| {
                if let Some(res) = squish.get(&(n.unit, n.cluster)) {
                    n.cluster = *res;
                }
            });
    }
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
    fn assemble_draft_graph(&self, c: &AssembleConfig) -> Graph {
        let (records, summaries) = assemble_draft(self, c);
        let nodes: Vec<_> = summaries
            .iter()
            .map(|s| {
                let id = s.id.clone();
                let segments: Vec<_> = s
                    .summary
                    .iter()
                    .map(|n| Tile {
                        unit: n.unit,
                        cluster: n.cluster,
                        strand: n.strand,
                    })
                    .collect();
                Node { id, segments }
            })
            .collect();
        let edges: Vec<_> = records
            .iter()
            .filter_map(|record| match &record.content {
                gfa::Content::Edge(e) => Some(e),
                _ => None,
            })
            .map(|e| {
                let from = e.sid1.id.to_string();
                let from_tail = e.sid1.is_forward();
                let to = e.sid2.id.to_string();
                let to_tail = e.sid2.is_forward();
                Edge {
                    from,
                    from_tail,
                    to,
                    to_tail,
                }
            })
            .collect();
        // for summary in summaries.iter() {
        //     debug!("SUMMARY\t{}", summary);
        // }
        Graph { nodes, edges }
    }
}

fn align_encoded_reads(ds: &DataSet, summaries: &[ContigSummary]) -> Vec<Vec<usize>> {
    // for (i, summary) in summaries.iter().enumerate() {
    //     let line: Vec<_> = summary
    //         .summary
    //         .iter()
    //         .map(|n| format!("{}-{}", n.unit, n.cluster))
    //         .collect();
    //     debug!("ALN\tCTG\t{i}\t{}", line.join("\t"));
    // }
    // let nodes: Vec<HashSet<_>> = summaries
    //     .iter()
    //     .map(|smy| smy.summary.iter().map(|n| (n.unit, n.cluster)).collect())
    //     .collect();
    // for read in ds.encoded_reads.iter() {
    //     let line: Vec<_> = read
    //         .nodes
    //         .iter()
    //         .map(|n| format!("{}-{}", n.unit, n.cluster))
    //         .collect();
    //     debug!("ALN\tREAD\t{}\t{}", read.id, line.join("\t"));
    //     let mut lightread: Vec<_> = read
    //         .nodes
    //         .iter()
    //         .map(|n| (n.unit, n.cluster, n.position_from_start, n.query_length()))
    //         .collect();
    //     while let Some((tid, start, end)) =
    //         search_aligned_region(&mut lightread, &nodes, &summaries)
    //     {
    //         debug!("ALN\tSEP\t{tid}\t{start}-{end}");
    //     }
    // }
    let nodes: Vec<HashSet<_>> = summaries
        .iter()
        .map(|smy| smy.summary.iter().map(|n| (n.unit, n.cluster)).collect())
        .collect();
    ds.encoded_reads
        .iter()
        .map(|read| {
            // let id = read.id;
            let read: Vec<_> = read.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
            let dist = distribute(&read, &nodes);
            // let line: Vec<_> = read.iter().map(|(n, c)| format!("{n}-{c}")).collect();
            // debug!("ALN\t{id}\t{}", line.join("\t"));
            // let line: Vec<_> = dist.iter().map(|d| format!("{d}")).collect();
            // debug!("ALN\t{id}\t{}", line.join("\t"));
            dist
        })
        .collect()
}

#[allow(dead_code)]
fn distribute(read: &[(u64, u64)], contigs: &[HashSet<(u64, u64)>]) -> Vec<usize> {
    let upperbound = read.len() as u64 + 10;
    // i -> j -> minimum number of switches from [0..i] position, ending with j-th contig.
    let mut dp = vec![vec![upperbound; contigs.len()]; read.len() + 1];
    dp[0].iter_mut().for_each(|x| *x = 0);
    // Traceback
    let mut tb = vec![vec![0; contigs.len()]; read.len() + 1];
    for (i, node) in read.iter().enumerate().map(|(i, x)| (i + 1, x)) {
        for (j, contig) in contigs.iter().enumerate() {
            let (argmin, min) = dp[i - 1]
                .iter()
                .enumerate()
                .map(|(idx, pen)| (idx, pen + (idx != j) as u64 + !contig.contains(node) as u64))
                .min_by_key(|x| x.1)
                .unwrap();
            dp[i][j] = min;
            tb[i][j] = argmin;
        }
    }
    let (mut argmin, _) = dp
        .last()
        .unwrap()
        .iter()
        .enumerate()
        .min_by_key(|x| x.1)
        .unwrap();
    let mut read_pos = read.len();
    let mut dists = vec![];
    while 0 < read_pos {
        dists.push(argmin);
        argmin = tb[read_pos][argmin];
        read_pos -= 1;
    }
    dists.reverse();
    assert_eq!(dists.len(), read.len());
    dists
}

/// ASSEMBLEIMPL
pub fn assemble(ds: &DataSet, c: &AssembleConfig) -> (Vec<gfa::Record>, Vec<ContigSummary>) {
    assert!(c.to_resolve);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let cov = ds.coverage.unwrap_or_else(|| panic!("Need coverage!"));
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, c);
    debug!("GRAPH\t{graph}");
    match ds.read_type {
        ReadType::CCS => graph.remove_lightweight_edges(1, true),
        ReadType::ONT | ReadType::None | ReadType::CLR => graph.remove_lightweight_edges(2, true),
    };
    graph.clean_up_graph_for_assemble(cov, &lens, &reads, c, ds.read_type);
    let (mut segments, mut edges, _, summaries) = graph.spell(c);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("{} segments({} bp in total).", segments.len(), total_base);
    if c.to_polish {
        let hmm = crate::model_tune::get_model(ds).unwrap();
        polish_segments(&mut segments, ds, &summaries, c, &ds.read_type, &hmm);
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
    // TODO: maybe just zip up segments and summaries would be OK?
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
pub fn assemble_draft(ds: &DataSet, c: &AssembleConfig) -> (Vec<gfa::Record>, Vec<ContigSummary>) {
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, c);
    graph.remove_lightweight_edges(2, true);
    let (segments, edge, group, summaries) = graph.spell(c);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("{} segments({} bp in total).", segments.len(), total_base);
    let nodes = segments.into_iter().map(|node| {
        let tags = match summaries.iter().find(|x| x.id == node.sid) {
            Some(contigsummary) => {
                let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
                let coverage =
                    gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
                vec![coverage]
            }
            None => Vec::new(),
        };
        gfa::Record::from_contents(gfa::Content::Seg(node), tags.into())
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags.into()));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![].into());
    let records: Vec<_> = std::iter::once(group).chain(nodes).chain(edges).collect();
    (records, summaries)
}

fn polish_segments(
    segments: &mut [gfa::Segment],
    ds: &DataSet,
    summaries: &[ContigSummary],
    c: &AssembleConfig,
    read_type: &definitions::ReadType,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
) {
    // Record/Associate reads to each segments.
    let raw_reads: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let mut fragments: Vec<Vec<&[u8]>> = vec![vec![]; summaries.len()];
    let dists = align_encoded_reads(ds, summaries);
    for (read, dist) in ds.encoded_reads.iter().zip(dists.iter()) {
        let seq = raw_reads[&read.id];
        if dist.len() == 0 {
            continue;
        }
        let mut pos = 0;
        let mut contig = dist[0];
        for (&d, n) in dist[1..].iter().zip(read.nodes[1..].iter()) {
            if d != contig {
                fragments[contig].push(&seq[pos..n.position_from_start]);
                pos = n.position_from_start + n.query_length();
                contig = d;
            }
        }
        fragments[contig].push(&seq[pos..]);
    }
    segments
        .iter_mut()
        .zip(fragments.iter())
        .for_each(|(seg, frag)| {
            polish_segment(seg, frag, c, read_type, hmm);
            polish_segment(seg, frag, c, read_type, hmm);
        });
}

// fn search_aligned_region(
//     read: &mut Vec<(u64, u64, usize, usize)>,
//     nodes: &[HashSet<(u64, u64)>],
//     summaries: &[ContigSummary],
// ) -> Option<(usize, usize, usize)> {
//     if read.is_empty() {
//         return None;
//     }
//     let (idx, _) = nodes
//         .iter()
//         .map(|ns| {
//             read.iter()
//                 .filter(|&&(u, c, _, _)| ns.contains(&(u, c)))
//                 .count()
//         })
//         .enumerate()
//         .max_by_key(|x| x.1)
//         .unwrap();
//     let (start_idx, end_idx) = align(read, &summaries[idx], &nodes[idx])?;
//     let start_pos = read[start_idx].2;
//     let end_pos = read[end_idx].2 + read[end_idx].3;
//     if read.len() - end_idx < start_idx {
//         *read = read.iter().take(start_idx).copied().collect();
//     } else {
//         *read = read.iter().skip(end_idx).copied().collect();
//     };
//     Some((idx, start_pos, end_pos))
// }

// fn align(
//     read: &[(u64, u64, usize, usize)],
//     summary: &ContigSummary,
//     scan: &HashSet<(u64, u64)>,
// ) -> Option<(usize, usize)> {
//     const OFFSET: usize = 4;
//     let ref_start = read
//         .iter()
//         .find(|&&(u, c, _, _)| scan.contains(&(u, c)))
//         .and_then(|&(u, c, _, _)| {
//             summary
//                 .summary
//                 .iter()
//                 .position(|sm| sm.unit == u && sm.cluster == c)
//         })
//         .unwrap();
//     let ref_end = read
//         .iter()
//         .rev()
//         .find(|&&(u, c, _, _)| scan.contains(&(u, c)))
//         .and_then(|&(u, c, _, _)| {
//             summary
//                 .summary
//                 .iter()
//                 .position(|sm| sm.unit == u && sm.cluster == c)
//         })
//         .unwrap();
//     let refr: Vec<_> = if ref_start < ref_end {
//         let ref_start = ref_start.saturating_sub(OFFSET);
//         let ref_end = (ref_end + OFFSET).min(summary.summary.len());
//         summary.summary[ref_start..ref_end]
//             .iter()
//             .map(|x| (x.unit, x.cluster))
//             .collect()
//     } else {
//         let ref_start = ref_end.saturating_sub(OFFSET);
//         let ref_end = (ref_start + OFFSET).min(summary.summary.len());
//         summary.summary[ref_start..ref_end]
//             .iter()
//             .map(|x| (x.unit, x.cluster))
//             .collect()
//     };
//     // Alignment.
//     // It is "infix-overlapping" alignment
//     let mut dp = vec![vec![0; refr.len() + 1]; read.len() + 1];
//     let (mut max, mut argmax) = (-1, (0, 0));
//     for (i, &(u, c, _, _)) in read.iter().enumerate() {
//         let q = (u, c);
//         let i = i + 1;
//         for (j, &r) in refr.iter().enumerate() {
//             let j = j + 1;
//             let mat_score = 2 * (r == q) as i32 - 1;
//             dp[i][j] = (dp[i - 1][j - 1] + mat_score)
//                 .max(dp[i - 1][j] - 1)
//                 .max(dp[i][j - 1] - 1)
//                 .max(0);
//             if max < dp[i][j] {
//                 (max, argmax) = (dp[i][j], (i, j));
//             }
//         }
//     }
//     // Traceback.
//     if max <= 1 {
//         return None;
//     }
//     let (mut qpos, mut rpos) = argmax;
//     let qend = qpos;
//     while dp[qpos][rpos] > 0 {
//         let r = refr[rpos - 1];
//         let q = read[qpos - 1];
//         let mat_score = 2 * (r == (q.0, q.1)) as i32 - 1;
//         let current = dp[qpos][rpos];
//         if current == dp[qpos - 1][rpos - 1] + mat_score {
//             qpos -= 1;
//             rpos -= 1;
//         } else if current == dp[qpos][rpos - 1] - 1 {
//             rpos -= 1;
//         } else if current == dp[qpos - 1][rpos - 1] - 1 {
//             qpos -= 1;
//         } else {
//             panic!();
//         }
//     }
//     Some((qpos, qend - 1))
// }

fn polish_segment(
    segment: &mut gfa::Segment,
    seqs: &[&[u8]],
    c: &AssembleConfig,
    rt: &definitions::ReadType,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
) {
    debug!("Aligning {} reads", seqs.len());
    let alns = match align_reads(segment, &seqs, rt, c) {
        Ok(res) => res,
        Err(why) => panic!("{:?}", why),
    };
    let seq = String::from_utf8(polish_by_chunking(&alns, segment, &seqs, c, rt, hmm)).unwrap();
    segment.slen = seq.len() as u64;
    segment.sequence = Some(seq);
}

fn align_reads(
    segment: &gfa::Segment,
    reads: &[&[u8]],
    read_type: &definitions::ReadType,
    c: &AssembleConfig,
) -> std::io::Result<kiley::sam::Sam> {
    let id = &segment.sid;
    let mut c_dir = std::env::current_dir()?;
    c_dir.push(format!("{}_polish", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir_all(&c_dir)?;
    // Create reference and reads.
    let (reference, reads) = {
        let mut reference = c_dir.clone();
        reference.push("segment.fa");
        use std::io::{BufWriter, Write};
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        writeln!(&mut wtr, ">{id}\n{}", segment.sequence.as_ref().unwrap())?;
        let mut reads_dir = c_dir.clone();
        reads_dir.push("reads.fa");
        let mut wtr = std::fs::File::create(&reads_dir).map(BufWriter::new)?;
        for (i, seq) in reads.iter().enumerate() {
            writeln!(&mut wtr, ">{}\n{}", i, std::str::from_utf8(seq).unwrap())?;
        }
        let reference = reference.into_os_string().into_string().unwrap();
        let reads = reads_dir.into_os_string().into_string().unwrap();
        (reference, reads)
    };
    let thr = format!("{}", c.threads);
    let mut args = vec!["-a", "-t", &thr, "--secondary=no"];
    match read_type {
        ReadType::CCS => args.extend(["-x", "map-hifi"]),
        ReadType::CLR => args.extend(["-x", "map-pb"]),
        ReadType::ONT => args.extend(["-x", "map-ont"]),
        ReadType::None => args.extend(["-x", "map-hifi"]),
    }
    let alignment = crate::minimap2::minimap2_args(&reference, &reads, &args);
    let alignment = kiley::sam::Sam::from_reader(std::io::BufReader::new(alignment.as_slice()));
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    Ok(alignment)
}

pub fn polish_by_chunking(
    alignments: &kiley::sam::Sam,
    segment: &gfa::Segment,
    reads: &[&[u8]],
    c: &AssembleConfig,
    read_type: &definitions::ReadType,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
) -> Vec<u8> {
    let template = match segment.sequence.as_ref() {
        Some(res) => kiley::SeqRecord::new(segment.sid.as_str(), res.as_bytes()),
        None => panic!(),
    };
    let reads: HashMap<_, _> = reads
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            let id = format!("{i}");
            let seq = kiley::SeqRecord::new(id, *seq);
            (format!("{i}"), seq)
        })
        .collect();
    let alignments: Vec<_> = alignments
        .records
        .iter()
        .filter_map(|aln| {
            let seq = reads.get(aln.q_name())?;
            Some((aln, seq))
        })
        .collect();
    let chunk_size = c.window_size;
    let radius = read_type.band_width(chunk_size);
    let max_cov = 50;
    let overlap = chunk_size / 5;
    let seed = 1071 * reads.len() as u64; // Arbitrary seed.
    let config =
        kiley::PolishConfig::with_model(radius, chunk_size, max_cov, overlap, seed, hmm.clone());
    kiley::polish_single(&template, &alignments, &config).seq
}

pub fn get_contig_copy_numbers(summaries: &[ContigSummary]) -> Vec<usize> {
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

// TODO:there are some bugs.
/// Return the co-occurence of the reads.
/// [i][j] -> # of reads shared by summaries[i] and summaries[j].
pub fn count_contig_connection(ds: &DataSet, summaries: &[ContigSummary]) -> Vec<Vec<u32>> {
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

#[cfg(test)]
pub mod tests {
    use super::*;
    #[test]
    fn distribute_test() {
        let read = vec![(0, 0), (1, 0), (2, 1), (3, 0)];
        let contig: HashSet<_> = vec![(0, 0), (1, 0), (2, 0), (3, 0)].into_iter().collect();
        let mut contigs = vec![contig];
        let dist = distribute(&read, &contigs);
        assert_eq!(dist, vec![0; 4]);
        let read = vec![(0, 0), (1, 1), (2, 1), (3, 1)];
        let contig: HashSet<_> = vec![(0, 1), (1, 1), (2, 1), (3, 1)].into_iter().collect();
        contigs.push(contig);
        let dist = distribute(&read, &contigs);
        assert_eq!(&dist[1..], &vec![1; 3]);
        let read = vec![(2, 1), (3, 1), (4, 1), (5, 1)];
        let contig: HashSet<_> = vec![(4, 1), (5, 1)].into_iter().collect();
        contigs.push(contig);
        let dist = distribute(&read, &contigs);
        assert_eq!(dist, vec![vec![1; 2], vec![2; 2]].concat());
    }
}

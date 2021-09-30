pub mod copy_number;
pub mod ditch_graph;
// pub mod string_graph;
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
    window_size: usize,
    // If true, resolve repeats.
    to_resolve: bool,
    min_span_reads: usize,
}

impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            to_polish: false,
            window_size: 100,
            to_resolve: false,
            min_span_reads: 6,
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
    ) -> Self {
        Self {
            window_size,
            threads,
            to_polish,
            to_resolve,
            min_span_reads,
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
}

impl Assemble for DataSet {
    fn squish_small_contig(&mut self, c: &AssembleConfig, len: usize) {
        {
            let reads: Vec<_> = self.encoded_reads.iter().collect();
            let cov = self.coverage.unwrap();
            let lens: Vec<_> = self.raw_reads.iter().map(|x| x.seq().len()).collect();
            let mut graph = DitchGraph::new(&reads, Some(&self.selected_chunks), c);
            graph.remove_lightweight_edges(2);
            graph.remove_zero_copy_elements(cov, &lens, 0.5);
            graph.resolve_repeats(&reads, c, 10f64);
            graph.assign_copy_number(cov, &lens);
            graph.resolve_repeats(&reads, c, 4f64);
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
    }
    fn assemble(&self, c: &AssembleConfig) -> GFA {
        debug!("Start assembly");
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let (records, summaries) = assemble(self, 0, c);
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
            debug!("{}\t{}\t{}", summary.id, copy_num, ids.join("\t"));
        }
        let mut header = vec![header];
        header.extend(records);
        GFA::from_records(header)
    }
    fn assemble_draft_graph(&self, c: &AssembleConfig) -> Graph {
        let mut cluster_and_num: HashMap<_, u32> = HashMap::new();
        for asn in self.assignments.iter() {
            *cluster_and_num.entry(asn.cluster).or_default() += 1;
        }
        debug!("There is {} clusters.", cluster_and_num.len());
        let (nodes, edges) = cluster_and_num
            .into_iter()
            .filter(|&(cl, num)| {
                if num < 10 {
                    debug!("Detected small group:{}(cluster:{})", num, cl);
                    false
                } else {
                    true
                }
            })
            .fold((vec![], vec![]), |(mut nodes, mut edges), (cl, _)| {
                let (records, summaries) = assemble(self, cl, c);
                nodes.extend(summaries.iter().map(|s| {
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
                }));
                edges.extend(
                    records
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
                        }),
                );
                for summary in summaries.iter() {
                    debug!("SUMMARY\t{}", summary);
                }
                (nodes, edges)
            });
        Graph { nodes, edges }
    }
}

pub fn assemble(
    ds: &DataSet,
    cl: usize,
    c: &AssembleConfig,
) -> (Vec<gfa::Record>, Vec<ContigSummary>) {
    let clusters: HashSet<_> = ds
        .assignments
        .iter()
        .filter(|asn| asn.cluster == cl)
        .map(|asn| asn.id)
        .collect();
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| clusters.contains(&r.id))
        .collect();
    if reads.is_empty() {
        panic!("Read is empty! {}", cl);
    }
    debug!("Constructing the {}-th ditch graph", cl);
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), c);
    // POINTER: ASSEMBLEIMPL
    // TODO: Tune this.
    graph.remove_lightweight_edges(2);
    if let Some(cov) = ds.coverage {
        if c.to_resolve {
            let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
            graph.remove_zero_copy_elements(cov, &lens, 0.5);
            graph.resolve_repeats(&reads, c, 10f64);
            graph.assign_copy_number(cov, &lens);
            graph.resolve_repeats(&reads, c, 4f64);

            // graph.zip_up_overclustering();
            // graph.remove_tips(0.5, 5);
            // graph.z_edge_selection();

            graph.assign_copy_number(cov, &lens);
        }
    }
    let (segments, edge, group, summaries) = graph.spell(c, cl);
    let total_base = segments.iter().map(|x| x.slen).sum::<u64>();
    debug!("{} segments({} bp in total).", segments.len(), total_base);
    let segments = if c.to_polish {
        segments
            .iter()
            .map(|segment| {
                let summary = summaries.iter().find(|s| s.id == segment.sid).unwrap();
                let segment = polish_segment(ds, segment, summary, c);
                polish_segment(ds, &segment, summary, c)
            })
            .collect()
    } else {
        segments
    };
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
    let records: Vec<_> = std::iter::once(group).chain(nodes).chain(edges).collect();
    (records, summaries)
}

fn polish_segment(
    ds: &DataSet,
    segment: &gfa::Segment,
    summary: &ContigSummary,
    c: &AssembleConfig,
) -> gfa::Segment {
    let reads = get_reads_in_cluster(ds, summary);
    debug!("Aligning {} reads", reads.len());
    let alignments = match align_reads(segment, &reads, c) {
        Ok(res) => res,
        Err(why) => panic!("{:?}", why),
    };
    let seq = String::from_utf8(polish_by_chunking(&alignments, segment, &reads, c)).unwrap();
    gfa::Segment::from(segment.sid.clone(), seq.len(), Some(seq))
}

fn get_reads_in_cluster<'a>(ds: &'a DataSet, summary: &ContigSummary) -> Vec<&'a [u8]> {
    let contained_unit: HashSet<(u64, u64)> = summary
        .summary
        .iter()
        .map(|elm| (elm.unit, elm.cluster))
        .collect();
    let ranges: HashMap<_, _> = ds
        .encoded_reads
        .iter()
        .filter_map(|r| {
            let mut nodes = r
                .nodes
                .iter()
                .filter(|n| contained_unit.contains(&(n.unit, n.cluster)));
            let start = nodes.next()?;
            let end = nodes.last()?;
            let start = start.position_from_start;
            let end = end.position_from_start + end.query_length();
            Some((r.id, (start, end)))
        })
        .collect();
    ds.raw_reads
        .iter()
        .filter_map(|r| ranges.get(&r.id).map(|&(s, e)| &r.seq()[s..e]))
        .collect()
}

fn align_reads(
    segment: &gfa::Segment,
    reads: &[&[u8]],
    c: &AssembleConfig,
) -> std::io::Result<kiley::sam::Sam> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut c_dir = std::env::current_dir()?;
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir)?;
    // Create reference and reads.
    let (reference, reads) = {
        let mut reference = c_dir.clone();
        reference.push("segment.fa");
        use std::io::{BufWriter, Write};
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        let seq = segment.sequence.as_ref().unwrap().as_bytes().to_vec();
        let seq = String::from_utf8_lossy(&seq);
        writeln!(&mut wtr, ">{}\n{}", &segment.sid, seq)?;
        let mut reads_dir = c_dir.clone();
        reads_dir.push("reads.fa");
        let mut wtr = std::fs::File::create(&reads_dir).map(BufWriter::new)?;
        for (i, seq) in reads.iter().enumerate() {
            writeln!(&mut wtr, ">{}\n{}", i, String::from_utf8_lossy(seq))?;
        }
        let reference = reference.into_os_string().into_string().unwrap();
        let reads = reads_dir.into_os_string().into_string().unwrap();
        (reference, reads)
    };
    let thr = format!("{}", c.threads);
    let args = [
        "-x",
        "map-pb",
        "-a",
        "-t",
        &thr,
        "--secondary=no",
        "-z",
        "600,400",
    ];
    let alignment = crate::minimap2::minimap2_args(&reference, &reads, &args);
    let alignment = kiley::sam::Sam::from_reader(std::io::BufReader::new(alignment.as_slice()));
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    debug!("Alignment done");
    Ok(alignment)
}

pub fn polish_by_chunking(
    alignments: &kiley::sam::Sam,
    segment: &gfa::Segment,
    reads: &[&[u8]],
    c: &AssembleConfig,
) -> Vec<u8> {
    let template_seq = match segment.sequence.as_ref() {
        Some(res) => res.as_bytes().to_vec(),
        None => panic!(),
    };
    let segment = (segment.sid.clone(), template_seq);
    debug!("Recording {} alignments...", alignments.records.len());
    let reads: Vec<_> = reads
        .iter()
        .enumerate()
        .map(|(i, r)| (format!("{}", i), r.to_vec()))
        .collect();
    use kiley::gphmm::*;
    let model = GPHMM::<Cond>::clr();
    let config = kiley::PolishConfig::with_model(100, c.window_size, 30, 50, 15, model);
    let mut polished = kiley::polish(&[segment], &reads, &alignments.records, &config);
    assert_eq!(polished.len(), 1);
    polished.pop().unwrap().1
}

type CpEdge = ((u64, u64), (u64, u64));
/// Estimate the copy number on each edge and each node.
pub fn copy_number_estimation(
    ds: &DataSet,
    c: &AssembleConfig,
) -> (HashMap<(u64, u64), usize>, HashMap<CpEdge, usize>) {
    let mut cluster_and_num: HashMap<_, u32> = HashMap::new();
    for asn in ds.assignments.iter() {
        *cluster_and_num.entry(asn.cluster).or_default() += 1;
    }
    let (mut node_cp, mut edge_cp) = (HashMap::new(), HashMap::new());
    for (cl, _) in cluster_and_num.into_iter() {
        let clusters: HashSet<_> = ds
            .assignments
            .iter()
            .filter(|asn| asn.cluster == cl)
            .map(|asn| asn.id)
            .collect();
        let reads: Vec<_> = ds
            .encoded_reads
            .iter()
            .filter(|r| clusters.contains(&r.id))
            .collect();
        let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), c);
        graph.remove_lightweight_edges(1);
        let cov = match ds.coverage {
            Some(res) => res,
            None => panic!(".coverage was checked before estimate coverage."),
        };
        let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
        let (node, edge) = graph.copy_number_estimation(cov, &lens);
        node_cp.extend(node);
        edge_cp.extend(edge.into_iter().map(|(((f, _), (t, _)), cp)| ((f, t), cp)));
    }
    (node_cp, edge_cp)
}

// let (_, edge_copy) = haplotyper::assemble::copy_number_estimation(&ds, &config);
// log::debug!("{} ZCE", edge_copy.values().filter(|&&x| x == 0).count());
// for read in ds.encoded_reads.iter_mut() {
//     let edge_to_cp = |w: &[Node]| {
//         let edge = if w[0].unit <= w[1].unit {
//             ((w[0].unit, w[0].cluster), (w[1].unit, w[1].cluster))
//         } else {
//             ((w[1].unit, w[1].cluster), (w[0].unit, w[0].cluster))
//         };
//         edge_copy[&edge]
//     };
//     let mut current = 0;
//     while let Some((start, _)) = read
//         .nodes
//         .windows(2)
//         .enumerate()
//         .skip(current)
//         .find(|(_, w)| edge_to_cp(w) == 0)
//     {
//         let (end, _) = read
//             .nodes
//             .windows(2)
//             .enumerate()
//             .skip(start)
//             .take_while(|(_, w)| edge_to_cp(w) == 0)
//             .last()
//             .unwrap();
//         if start < end {
//             let removed: Vec<_> = read.nodes[start..=end]
//                 .iter()
//                 .map(|x| format!("{}", x.unit))
//                 .collect();
//             log::debug!("Removing {:?}", removed);
//             for _ in start..=end {
//                 read.remove(start);
//             }
//         }
//         current = start + 1;
//     }
// }

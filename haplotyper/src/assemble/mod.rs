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
    /// Zip up suspicious haplotig. In other words,
    /// after usual assembly workflow, if two haplotig shares reads more than `count`,
    /// it indicates that these contigs can not be phased because of something. We squish such haplotigs with length less than `len`.
    fn zip_up_suspicious_haplotig(&mut self, c: &AssembleConfig, count: u32, len: usize);
}

impl Assemble for DataSet {
    fn zip_up_suspicious_haplotig(&mut self, c: &AssembleConfig, count: u32, len: usize) {
        let (_, summaries) = assemble(self, 0, c);
        // Select unique units.
        // let unique_units = filter_non_unique_units(&summaries);
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
        let mut graph = DitchGraph::new(&reads, Some(&self.selected_chunks), c);
        graph.remove_lightweight_edges(2, true);
        graph.clean_up_graph_for_assemble(cov, &lens, &reads, c);
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
        let reads: Vec<_> = self.encoded_reads.iter().collect();
        let cov = self.coverage.unwrap();
        let lens: Vec<_> = self.raw_reads.iter().map(|x| x.seq().len()).collect();
        let mut graph = DitchGraph::new(&reads, Some(&self.selected_chunks), c);
        graph.remove_lightweight_edges(2, true);
        graph.assign_copy_number(cov, &lens);
        graph.remove_zero_copy_elements(&lens, 0.5);
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
        let (records, summaries) = assemble(self, 0, c);
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
        let header = gfa::Record::from_contents(header, vec![]);
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
    graph.remove_lightweight_edges(2, true);
    // let cov = ds.coverage.unwrap_or_else(|| panic!("Need coverage!"));
    //     if c.to_resolve {
    //         let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    //         graph.clean_up_graph_for_assemble(cov, &lens, &reads, c);
    //     }
    graph.remove_lightweight_edges(1, true);
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

// type CpEdge = ((u64, u64), (u64, u64));
// /// Estimate the copy number on each edge and each node.
// pub fn copy_number_estimation(
//     ds: &DataSet,
//     c: &AssembleConfig,
// ) -> (HashMap<(u64, u64), usize>, HashMap<CpEdge, usize>) {
//     let mut cluster_and_num: HashMap<_, u32> = HashMap::new();
//     for asn in ds.assignments.iter() {
//         *cluster_and_num.entry(asn.cluster).or_default() += 1;
//     }
//     let (mut node_cp, mut edge_cp) = (HashMap::new(), HashMap::new());
//     for (cl, _) in cluster_and_num.into_iter() {
//         let clusters: HashSet<_> = ds
//             .assignments
//             .iter()
//             .filter(|asn| asn.cluster == cl)
//             .map(|asn| asn.id)
//             .collect();
//         let reads: Vec<_> = ds
//             .encoded_reads
//             .iter()
//             .filter(|r| clusters.contains(&r.id))
//             .collect();
//         let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), c);
//         graph.remove_lightweight_edges(2, true);
//         let cov = match ds.coverage {
//             Some(res) => res,
//             None => panic!("coverage was checked before estimate coverage."),
//         };
//         let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
//         let (node, edge) = graph.copy_number_estimation(cov, &lens);
//         node_cp.extend(node);
//         edge_cp.extend(edge.into_iter().map(|(((f, _), (t, _)), cp)| ((f, t), cp)));
//     }
//     (node_cp, edge_cp)
// }

// pub fn filter_non_unique_units(summaries: &[ContigSummary]) -> Vec<Vec<(u64, u64)>> {
//     let mut unit_counts: HashMap<_, u32> = HashMap::new();
//     for summary in summaries.iter() {
//         for node in summary.summary.iter() {
//             *unit_counts.entry((node.unit, node.cluster)).or_default() += 1;
//         }
//     }
//     // Retain only unique units.
//     unit_counts.retain(|_, val| *val == 1);
//     summaries
//         .iter()
//         .map(|summary| {
//             summary
//                 .summary
//                 .iter()
//                 .map(|n| (n.unit, n.cluster))
//                 .filter(|key| unit_counts.contains_key(&key))
//                 .collect()
//         })
//         .collect()
// }

pub fn get_contig_copy_numbers(summaries: &[ContigSummary]) -> Vec<usize> {
    summaries
        .iter()
        .map(|summary| {
            let (cp, num) = summary
                .summary
                .iter()
                .filter_map(|n| n.copy_number)
                .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
            (cp as f64 / num as f64).round() as usize
        })
        .collect()
}

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
    for ((u, c), xs) in contig_terminals.iter() {
        for &i in xs.iter() {
            trace!("TERMINALS\t{}\t{}\t{}", summaries[i].id, u, c);
        }
    }
    let mut shared_reads = vec![vec![0; summaries.len()]; summaries.len()];
    for read in ds.encoded_reads.iter() {
        let mut read_through = vec![];
        for (i, w) in read.nodes.windows(2).enumerate() {
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
                if i == 0 {
                    read_through.push(f);
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

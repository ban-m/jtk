pub mod copy_number;
pub mod ditch_graph;
pub mod string_graph;
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
}

impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            to_polish: false,
            window_size: 100,
        }
    }
}
impl AssembleConfig {
    pub fn new(threads: usize, window_size: usize, to_polish: bool) -> Self {
        Self {
            window_size,
            threads,
            to_polish,
        }
    }
}

pub trait Assemble {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA;
    fn assemble_as_graph(&self, c: &AssembleConfig) -> Vec<Graph>;
}

impl Assemble for DataSet {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA {
        let mut cluster_and_num: HashMap<_, u32> = HashMap::new();
        for asn in self.assignments.iter() {
            *cluster_and_num.entry(asn.cluster).or_default() += 1;
        }
        debug!("There is {} clusters.", cluster_and_num.len());
        debug!("Start assembly");
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let mut cluster_and_num: Vec<_> = cluster_and_num.into_iter().collect();
        cluster_and_num.sort_by_key(|x| x.1);
        cluster_and_num.reverse();
        let records: Vec<_> = cluster_and_num
            .into_iter()
            .filter(|&(cl, num)| {
                if num < 10 {
                    debug!("Detected small group:{}(cluster:{})", num, cl);
                    false
                } else {
                    true
                }
            })
            .flat_map(|(cl, _)| {
                let (records, _summaries) = assemble(self, cl, c);
                records
            })
            .collect();
        let mut header = vec![header];
        header.extend(records);
        GFA::from_records(header)
    }
    fn assemble_as_graph(&self, c: &AssembleConfig) -> Vec<Graph> {
        let mut cluster_and_num: HashMap<_, u32> = HashMap::new();
        for asn in self.assignments.iter() {
            *cluster_and_num.entry(asn.cluster).or_default() += 1;
        }
        debug!("There is {} clusters.", cluster_and_num.len());
        cluster_and_num
            .into_iter()
            .filter(|&(cl, num)| {
                if num < 10 {
                    debug!("Detected small group:{}(cluster:{})", num, cl);
                    false
                } else {
                    true
                }
            })
            .map(|(cl, _)| {
                let (records, summaries) = assemble(&self, cl, c);
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
                let edges = records
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
                for summary in summaries.iter() {
                    debug!("SUMMARY\t{}", summary);
                }
                Graph { nodes, edges }
            })
            .collect()
    }
}

fn assemble(
    ds: &DataSet,
    cl: usize,
    c: &AssembleConfig,
) -> (
    Vec<gfa::Record>,
    // Vec<gfa::Segment>,
    // Vec<gfa::Edge>,
    // gfa::Group,
    Vec<ContigSummary>,
) {
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
    debug!("Constructing the {}-th ditch graph", cl);
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), c);
    let short_reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.len() < 4)
        .collect();
    graph.remove_edges_from_short_reads(&short_reads);
    // graph.resolve_repeats();
    // for i in 0..2 {
    //     graph.remove_lightweight_edges(i);
    // graph.collapse_bubble(c);
    // }
    // graph.remove_tips();
    // graph.remove_small_component(5);
    // graph.remove_lightweight_edges(3);
    // graph.collapse_bubble(c);
    // debug!("{}", graph);
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
    let nodes = segments.into_iter().map(|node| {
        let tags = summaries
            .iter()
            .find(|x| x.id == node.sid)
            .map(|contigsummary| {
                let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
                vec![gfa::SamTag {
                    inner: format!("cv:i:{}", total / contigsummary.summary.len()),
                }]
            })
            .unwrap_or(vec![]);
        gfa::Record::from_contents(gfa::Content::Seg(node), tags)
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    let records: Vec<_> = std::iter::once(group).chain(nodes).chain(edges).collect();
    (records, summaries)
    // (segments, edge, group, summaries)
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
    let seq = String::from_utf8(polish_by_chunking(&alignments, &segment, &reads, c)).unwrap();
    gfa::Segment::from(segment.sid.clone(), seq.len(), Some(seq))
}

fn get_reads_in_cluster<'a>(ds: &'a DataSet, summary: &ContigSummary) -> Vec<&'a RawRead> {
    let contained_unit: HashSet<(u64, u64)> = summary
        .summary
        .iter()
        .map(|elm| (elm.unit, elm.cluster))
        .collect();
    let ids: HashSet<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| {
            r.nodes
                .iter()
                .any(|n| contained_unit.contains(&(n.unit, n.cluster)))
        })
        .map(|r| r.id)
        .collect();
    ds.raw_reads
        .iter()
        .filter(|r| ids.contains(&r.id))
        .collect()
}

fn align_reads(
    segment: &gfa::Segment,
    reads: &[&RawRead],
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
        for read in reads.iter() {
            writeln!(&mut wtr, ">{}\n{}", read.name, read.seq)?;
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
    reads: &[&RawRead],
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
        .map(|r| (r.name.clone(), r.seq().to_vec()))
        .collect();
    use kiley::gphmm::*;
    let model = GPHMM::<Cond>::clr();
    let config = kiley::PolishConfig::with_model(100, c.window_size, 30, 50, 15, model);
    let mut polished = kiley::polish(&vec![segment], &reads, &alignments.records, &config);
    assert_eq!(polished.len(), 1);
    polished.pop().unwrap().1
}

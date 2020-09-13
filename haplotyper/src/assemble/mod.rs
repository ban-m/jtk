mod ditch_graph;
use definitions::*;
use ditch_graph::*;
use gfa::GFA;
use serde::*;
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

use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct AssembleConfig {}
impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {}
    }
}

pub trait Assemble {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA;
    fn assemble_as_graph(&self, c: &AssembleConfig) -> Vec<Graph>;
}

impl Assemble for DataSet {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA {
        let mut clusters: HashMap<_, Vec<_>> = HashMap::new();
        {
            // ID=>Cluster
            let id2cluster: HashMap<u64, usize> = self
                .assignments
                .iter()
                .map(|asn| (asn.id, asn.cluster))
                .collect();
            for read in self.encoded_reads.iter() {
                if let Some(&cluster) = id2cluster.get(&read.id) {
                    clusters.entry(cluster).or_default().push(read);
                }
            }
        }
        debug!("There is {} clusters.", clusters.len());
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        debug!("Start assembly");
        let records: Vec<_> = clusters
            .into_iter()
            .flat_map(|(cl, reads)| cluster_to_gfa(cl, reads, c))
            .chain(std::iter::once(header))
            .collect();
        GFA::from_records(records)
    }
    fn assemble_as_graph(&self, c: &AssembleConfig) -> Vec<Graph> {
        let mut clusters: HashMap<_, Vec<_>> = HashMap::new();
        {
            let id2cluster: HashMap<u64, usize> = self
                .assignments
                .iter()
                .map(|asn| (asn.id, asn.cluster))
                .collect();
            for read in self.encoded_reads.iter() {
                if let Some(&cluster) = id2cluster.get(&read.id) {
                    clusters.entry(cluster).or_default().push(read);
                }
            }
        }
        debug!("There is {} clusters.", clusters.len());
        debug!("Start assembly");
        clusters
            .into_iter()
            .filter_map(|(cl, reads)| cluster_to_graph(cl, reads, c))
            .collect()
    }
}

fn cluster_to_graph(cl: usize, reads: Vec<&EncodedRead>, c: &AssembleConfig) -> Option<Graph> {
    debug!("Constructing the {}-th ditch graph", cl);
    if reads.len() < 10 {
        debug!("Detected small group:{}", reads.len());
        debug!("The reads below are discarded.");
        debug!("ID:Length:Nodes");
        for read in reads.iter() {
            debug!("{}:{}:{}", read.id, read.original_length, read.nodes.len());
        }
        return None;
    }
    let mut graph = DitchGraph::new(&reads, c);
    graph.resolve_repeats();
    graph.remove_redundant_edges(2);
    graph.remove_tips();
    graph.collapse_buddle(c);
    graph.remove_small_component(5);
    let (_, edges, _, summaries) = graph.spell(c, cl);
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
    let edges = edges
        .iter()
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
        debug!("{}", summary);
    }
    let graph = Graph { nodes, edges };
    Some(graph)
}

fn cluster_to_gfa(cl: usize, reads: Vec<&EncodedRead>, c: &AssembleConfig) -> Vec<gfa::Record> {
    debug!("Constructing the {}-th ditch graph", cl);
    if reads.len() < 10 {
        debug!("Detected small group:{}", reads.len());
        debug!("The reads below are discarded.");
        debug!("ID:Length:Nodes");
        for read in reads.iter() {
            debug!("{}:{}:{}", read.id, read.original_length, read.nodes.len());
        }
        return vec![];
    }
    let mut graph = DitchGraph::new(&reads, c);
    graph.resolve_repeats();
    graph.remove_redundant_edges(2);
    graph.remove_tips();
    graph.collapse_buddle(c);
    graph.remove_small_component(5);
    debug!("{}", graph);
    let mut records = vec![];
    let (nodes, edges, group, summaries) = graph.spell(c, cl);
    for summary in summaries {
        debug!("{}", summary);
    }
    let nodes = nodes
        .into_iter()
        .map(gfa::Content::Seg)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(nodes);
    let edges = edges
        .into_iter()
        .map(gfa::Content::Edge)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(edges);
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    records.push(group);
    records
}

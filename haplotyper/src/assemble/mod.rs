use definitions::*;
use gfa::GFA;
use std::collections::HashMap;
mod ditch_graph;
use ditch_graph::*;
#[derive(Debug, Clone)]
pub struct AssembleConfig {}
impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {}
    }
}

pub trait Assemble {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA;
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
    debug!("{}", graph);
    let mut records = vec![];
    let (nodes, edges, group) = graph.spell(c, cl);
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

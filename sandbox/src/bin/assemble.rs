#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|n| n.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    let config = AssembleConfig::new(1, 100, false, true, 3, 4f64);
    let records = assemble_draft(&ds, &config);
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![]);
    let mut header = vec![header];
    header.extend(records);
    let gfa = gfa::GFA::from_records(header);
    println!("{}", gfa);
    Ok(())
}

pub fn assemble_draft(ds: &DataSet, c: &AssembleConfig) -> Vec<gfa::Record> {
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, c);
    // graph.remove_lightweight_edges(2, true);
    let (segments, edge, group, summaries) = graph.spell(c);
    let nodes = segments.into_iter().map(|node| {
        let tags = match summaries.iter().find(|x| x.id == node.sid) {
            Some(contigsummary) => {
                let ids: Vec<_> = contigsummary
                    .summary
                    .iter()
                    .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
                    .collect();
                let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
                let coverage =
                    gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
                log::debug!(
                    "CONUNIT\t{}\t{}\t{}",
                    contigsummary.id,
                    total / contigsummary.summary.len(),
                    ids.join("\t")
                );
                vec![coverage]
            }
            None => Vec::new(),
        };
        gfa::Record::from_contents(gfa::Content::Seg(node), tags)
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    std::iter::once(group).chain(nodes).chain(edges).collect()
}

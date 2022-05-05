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
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|n| n.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    let config = AssembleConfig::new(1, 100, true, true, 3, 4f64);
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
    for read in ds.encoded_reads.iter() {
        let len = read.recover_raw_read().len();
        assert_eq!(read.original_length, len);
    }
    eprintln!("{graph}");
    graph.remove_lightweight_edges(2, true);
    use log::*;
    for node in graph.nodes() {
        let (unit, cl) = node.node;
        let seq = std::str::from_utf8(node.seq()).unwrap();
        debug!("CONS\t>N{unit}-{cl}\nCONS\t{seq}");
        for edge in node.edges.iter().filter(|e| e.from <= e.to) {
            let seq = match edge.label() {
                Some(res) => std::str::from_utf8(res).unwrap(),
                None => continue,
            };
            let (funit, fcl) = edge.from;
            let (tunit, tcl) = edge.to;
            let id = format!("{funit}-{fcl}-{tunit}-{tcl}");
            debug!("CONS\t>E{id}\nCONS\t{seq}");
        }
    }
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

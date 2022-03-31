#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::encode::deletion_fill::*;
    // let (read_error_rate, unit_error_rate, deviation) = estimate_error_rate_dev(&ds, 0.35);
    let (read_error_rate, unit_error_rate, deviation) = estimate_error_rate(&ds, 0.35);
    eprintln!("{deviation}");
    for read in ds.encoded_reads.iter() {
        println!("READ\t{}\t{}", read.id, read_error_rate[read.id as usize]);
    }
    for unit in ds.selected_chunks.iter() {
        println!("UNIT\t{}\t{}", unit.id, unit_error_rate[unit.id as usize]);
        // for cl in 0..unit.cluster_num {
        //     let error = unit_error_rate[unit.id as usize][cl];
        //     println!("UNIT\t{}\t{cl}\t{error}", unit.id);
        // }
    }
    // let header = gfa::Content::Header(gfa::Header::default());
    // let header = gfa::Record::from_contents(header, vec![]);
    // let config = AssembleConfig::new(1, 100, false, true, 3, 4f64);
    // let records = assemble_draft(&ds, &config);
    // let records = std::iter::once(header.clone()).chain(records).collect();
    // let gfa = gfa::GFA::from_records(records);
    // use haplotyper::purge_diverged::*;
    // use std::io::{BufWriter, Write};
    // if let Ok(mut wtr) = std::fs::File::create("before.gfa").map(BufWriter::new) {
    //     writeln!(&mut wtr, "{}", gfa)?;
    // }
    // let purge_config = PurgeDivConfig::new(56);
    // {
    //     let mut ds = ds.clone();
    //     ds.purge(&purge_config);
    //     let records = assemble_draft(&ds, &config);
    //     let records = std::iter::once(header.clone()).chain(records).collect();
    //     let gfa = gfa::GFA::from_records(records);
    //     if let Ok(mut wtr) = std::fs::File::create("old.gfa").map(BufWriter::new) {
    //         writeln!(&mut wtr, "{}", gfa)?;
    //     }
    // }
    // {
    //     let mut ds = ds.clone();
    //     ds.purge_dev(&purge_config);
    //     let records = assemble_draft(&ds, &config);
    //     let records = std::iter::once(header).chain(records).collect();
    //     let gfa = gfa::GFA::from_records(records);
    //     if let Ok(mut wtr) = std::fs::File::create("new.gfa").map(BufWriter::new) {
    //         writeln!(&mut wtr, "{}", gfa)?;
    //     }
    // }
    Ok(())
}

// pub fn assemble_draft(ds: &DataSet, c: &AssembleConfig) -> Vec<gfa::Record> {
//     let reads: Vec<_> = ds.encoded_reads.iter().collect();
//     use haplotyper::assemble::ditch_graph::DitchGraph;
//     let graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, c);
//     // graph.remove_lightweight_edges(2, true);
//     let (segments, edge, group, summaries) = graph.spell(c);
//     let nodes = segments.into_iter().map(|node| {
//         let tags = match summaries.iter().find(|x| x.id == node.sid) {
//             Some(contigsummary) => {
//                 let ids: Vec<_> = contigsummary
//                     .summary
//                     .iter()
//                     .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
//                     .collect();
//                 let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
//                 let coverage =
//                     gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
//                 log::debug!(
//                     "CONUNIT\t{}\t{}\t{}",
//                     contigsummary.id,
//                     total / contigsummary.summary.len(),
//                     ids.join("\t")
//                 );
//                 vec![coverage]
//             }
//             None => Vec::new(),
//         };
//         gfa::Record::from_contents(gfa::Content::Seg(node), tags)
//     });
//     let edges = edge
//         .into_iter()
//         .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
//     let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
//     std::iter::once(group).chain(nodes).chain(edges).collect()
// }

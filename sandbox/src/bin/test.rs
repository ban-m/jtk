use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::*;
    let config = haplotyper::AssembleConfig::new(24, 2000, false);
    let coverage = ds.coverage.unwrap();
    let lens: Vec<_> = ds.raw_reads.iter().map(|read| read.seq().len()).collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    graph.remove_lightweight_edges(1);
    // First clean up.
    graph.remove_zero_copy_elements(coverage, &lens, 0.51);
    graph.z_edge_selection();
    graph.remove_tips(0.5, 5);
    graph.zip_up_overclustering();
    graph.resolve_repeats(&reads, &config, 8f64);
    // Second clean up.
    graph.assign_copy_number(coverage, &lens);
    graph.resolve_repeats(&reads, &config, 4f64);
    // Third clustering.
    graph.assign_copy_number(coverage, &lens);
    graph.z_edge_selection();
    graph.zip_up_overclustering();
    graph.resolve_repeats(&reads, &config, 4f64);
    // Assemble
    let gfa = assemble_as_gfa(&graph, &config);
    println!("{}", gfa);
    Ok(())
}

#[allow(dead_code)]
fn assemble_as_gfa(
    graph: &haplotyper::assemble::ditch_graph::DitchGraph,
    c: &haplotyper::AssembleConfig,
) -> gfa::GFA {
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![]);
    let (segments, edge, group, summaries) = graph.spell(c, 0);
    for summary in summaries.iter() {
        let ids: Vec<_> = summary
            .summary
            .iter()
            .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
            .collect();
        log::debug!("{}\t{}", summary.id, ids.join("\t"));
    }
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
            .unwrap_or(vec![]);
        gfa::Record::from_contents(gfa::Content::Seg(node), tags)
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    let records: Vec<_> = std::iter::once(group).chain(nodes).chain(edges).collect();
    let mut header = vec![header];
    header.extend(records);
    gfa::GFA::from_records(header)
}

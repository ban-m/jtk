use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    let coverage = ds.coverage.unwrap();
    let lens: Vec<_> = ds.raw_reads.iter().map(|read| read.seq().len()).collect();
    let config = haplotyper::AssembleConfig::new(24, 2000, false);
    let mut gfa = assemble_as_gfa(&ds, &config);
    haplotyper::assemble::copy_number::estimate_copy_number_on_gfa(&mut gfa, coverage, &lens, 2000);
    println!("{}", gfa);
    Ok(())
}

fn assemble_as_gfa(ds: &DataSet, c: &haplotyper::AssembleConfig) -> gfa::GFA {
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![]);
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), c);
    graph.remove_lightweight_edges(1);
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
    let mut header = vec![header];
    header.extend(records);
    gfa::GFA::from_records(header)
}

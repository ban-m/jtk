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
    let to_erase: bool = args.get(2).map(|x| x == "true").unwrap_or(false);
    if to_erase {
        ds.encoded_reads
            .iter_mut()
            .flat_map(|n| n.nodes.iter_mut())
            .for_each(|n| n.cluster = 0);
    }
    let config = AssembleConfig::new(1, 100, false, false, 3, 4f64);
    let records = assemble_draft(&ds, &config);
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![].into());
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
    graph.remove_lightweight_edges(3, true);
    let cov = ds.coverage.unwrap_or(30.0);
    let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(4395);
    // let lens: Vec<_> = ds.encoded_reads.iter().map(|r| r.original_length).collect();
    // graph.assign_copy_number(cov, &lens);
    graph.assign_copy_number_mst(cov, &mut rng);
    // use log::*;
    // for node in graph.nodes() {
    //     let (unit, cl) = node.node;
    //     let seq = std::str::from_utf8(node.seq()).unwrap();
    //     debug!("CONS\t>N{unit}-{cl}\nCONS\t{seq}");
    //     for edge in node.edges.iter().filter(|e| e.from <= e.to) {
    //         let seq = match edge.label() {
    //             Some(res) => std::str::from_utf8(res).unwrap(),
    //             None => continue,
    //         };
    //         let (funit, fcl) = edge.from;
    //         let (tunit, tcl) = edge.to;
    //         let id = format!("{funit}-{fcl}-{tunit}-{tcl}");
    //         debug!("CONS\t>E{id}\nCONS\t{seq}");
    //     }
    // }
    let (segments, edge, _, summaries) = graph.spell(c);
    let mut groups: HashMap<_, Vec<_>> = HashMap::new();
    let nodes: Vec<_> = segments
        .into_iter()
        .map(|node| {
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
                    let (cp, cpnum) = contigsummary
                        .summary
                        .iter()
                        .filter_map(|elm| elm.copy_number)
                        .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
                    let copy_number = cp / cpnum.max(1);
                    let ids = ids.join("\t");
                    log::debug!("CONUNIT\t{}\t{copy_number}\t{ids}", contigsummary.id,);
                    groups
                        .entry(copy_number)
                        .or_default()
                        .push(node.sid.clone());
                    let copy_number = gfa::SamTag::new(format!("cp:i:{copy_number}"));
                    vec![coverage, copy_number]
                }
                None => Vec::new(),
            };
            gfa::Record::from_contents(gfa::Content::Seg(node), tags.into())
        })
        .collect();
    let groups = groups.into_iter().map(|(cp, ids)| {
        let group = gfa::UnorderedGroup {
            uid: Some(format!("cp:i:{}", cp)),
            ids,
        };
        let group = gfa::Content::Group(gfa::Group::Set(group));
        gfa::Record::from_contents(group, vec![].into())
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags.into()));
    groups.chain(nodes).chain(edges).collect()
}

use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use std::collections::{HashMap, HashSet};
    let mut counts: HashMap<_, (usize, i64)> = HashMap::new();
    let normalize = |from: &Node, to: &Node| -> (u64, u64) {
        match from.unit <= to.unit {
            true => (from.unit, to.unit),
            false => (to.unit, from.unit),
        }
    };
    for read in ds.encoded_reads.iter() {
        for (i, from) in read.nodes.iter().enumerate().take(read.nodes.len() - 1) {
            let to = &read.nodes[i + 1];
            let entry = counts.entry(normalize(from, to)).or_default();
            entry.0 += 1;
            entry.1 += to.position_from_start as i64
                - (from.position_from_start - from.query_length()) as i64;
        }
    }
    // Convert to the "calibrated occurance".
    use haplotyper::assemble::copy_number::CoverageCalibrator;
    let lens: Vec<_> = ds.encoded_reads.iter().map(|r| r.original_length).collect();
    let calib = CoverageCalibrator::new(&lens);
    let calib_coverage: HashMap<_, f64> = counts
        .iter()
        .map(|(&(f, t), &(obs, totlen))| {
            let len = (totlen / obs as i64).max(0) as usize;
            let calibed = calib.calib(obs, len);
            eprintln!("DUMP\t{}\t{}\t{}\t{}\t{}", f, t, obs, len, calibed);
            ((f, t), calibed)
        })
        .collect();
    let median = {
        let mut covs: Vec<_> = calib_coverage.values().collect();
        let len = covs.len() / 2;
        let (_, median, _) = covs.select_nth_unstable_by(len, |x, y| x.partial_cmp(y).unwrap());
        *median
    };
    eprintln!("{}", median);
    // (edge)->Option<The element to be removed>.
    let to_remove: HashMap<_, u64> = calib_coverage
        .iter()
        .filter(|(_, &cov)| cov < median / 4f64)
        .filter_map(|(&key, &cov)| {
            // Search reads with (from,fd,to,td) occurence.
            let mut former_neighbor = HashSet::new();
            let mut later_neighbor = HashSet::new();
            for read in ds.encoded_reads.iter() {
                for (i, from) in read.nodes.iter().enumerate().take(read.nodes.len() - 1) {
                    let to = &read.nodes[i + 1];
                    if normalize(from, to) == key {
                        if let Some(next) = read.nodes.get(i + 2) {
                            match from.unit <= to.unit {
                                true => former_neighbor.insert(next.unit),
                                false => later_neighbor.insert(next.unit),
                            };
                        }
                        if let Some(prev) = read.nodes.get(i - 1) {
                            match from.unit <= to.unit {
                                true => later_neighbor.insert(prev.unit),
                                false => former_neighbor.insert(prev.unit),
                            };
                        }
                    }
                }
            }
            for next_unit in former_neighbor {
                let probe = if key.0 <= next_unit {
                    (key.0, next_unit)
                } else {
                    (next_unit, key.0)
                };
                let mod_cov = *calib_coverage.get(&probe).unwrap_or(&0f64);
                if key == (1240, 1479) {
                    eprintln!("II\t{}\t{}\t{}", next_unit, mod_cov, cov);
                }
                if 3f64 * cov < mod_cov {
                    return Some((key, key.1));
                }
            }
            for prev_unit in later_neighbor {
                let probe = if key.1 <= prev_unit {
                    (key.1, prev_unit)
                } else {
                    (prev_unit, key.1)
                };
                let mod_cov = *calib_coverage.get(&probe).unwrap_or(&0f64);
                if 3f64 * cov < mod_cov {
                    return Some((key, key.0));
                }
            }
            None
        })
        .collect();
    for ((f, t), rm) in to_remove.iter() {
        eprintln!("REMOVING\t{}\t{}\t{}", f, t, rm);
    }
    let mut ds = ds;
    for read in ds.encoded_reads.iter_mut() {
        let remove_idx: Vec<_> = read
            .nodes
            .iter()
            .enumerate()
            .take(read.nodes.len() - 1)
            .filter_map(|(i, from)| {
                let to = &read.nodes[i + 1];
                to_remove
                    .get(&normalize(from, to))
                    .map(|&removing| match removing == from.unit {
                        true => i,
                        false => i + 1,
                    })
            })
            .collect();
        let mut offset = 0;
        for i in remove_idx {
            read.remove(i - offset);
            offset += 1;
        }
    }
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // use std::collections::HashMap;
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // use haplotyper::assemble::re_clustering;
    // use haplotyper::em_correction::ClusteringCorrection;
    // let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    // let ds = ds.correct_clustering_em(5, 5, 3);
    // let ds = re_clustering(ds, &config);
    // let ds = ds.correct_clustering_em(5, 5, 3);
    // println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     let ans = id2desc[&read.id].contains("000252v2") as usize;
    //     for node in read.nodes.iter() {
    //         counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
    //     }
    // }
    // for ((unit, cluster), counts) in counts.iter() {
    //     println!("{}\t{}\t{}\t{}", unit, cluster, counts[0], counts[1]);
    // }
    // let tigs = std::fs::File::open(&args[2]).map(BufReader::new)?;
    // use std::io::BufRead;
    // let tigs: Vec<_> = tigs
    //     .lines()
    //     .filter_map(|x| x.ok())
    //     .map(|line| {
    //         let mut line = line.split('\t');
    //         let id = line.next().unwrap().to_string();
    //         let clusters: Vec<_> = line
    //             .filter_map(|x| {
    //                 let mut elm = x.split('-');
    //                 let unit: u64 = elm.next()?.parse().unwrap();
    //                 let cluster: u64 = elm.next()?.parse().unwrap();
    //                 Some((unit, cluster))
    //             })
    //             .collect();
    //         (id, clusters)
    //     })
    //     .collect();
    // let mut counts: HashMap<String, [u32; 2]> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     let is_in_hapa = id2desc[&read.id].contains("000252v2") as usize;
    //     if let Some((tig_id, _)) = tigs
    //         .iter()
    //         .filter_map(|(tig, clusters)| {
    //             let count = read
    //                 .nodes
    //                 .iter()
    //                 .filter(|n| clusters.contains(&(n.unit, n.cluster)))
    //                 .count();
    //             (count != 0).then(|| (tig, count))
    //         })
    //         .max_by_key(|x| x.0)
    //     {
    //         counts.entry(tig_id.clone()).or_default()[is_in_hapa] += 1;
    //     }
    // }
    // for (id, c) in counts {
    //     let length = tigs.iter().find(|x| x.0 == id).map(|x| x.1.len()).unwrap();
    //     println!("{}\t{}\t{}\t{}", id, length, c[0], c[1]);
    // }
    // ds.assignments = ds
    //     .encoded_reads
    //     .iter()
    //     .map(|r| Assignment::new(r.id, 0))
    //     .collect();
    // let reads: Vec<_> = ds.encoded_reads.iter().collect();
    // use haplotyper::assemble::ditch_graph::*;
    // let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    // let coverage = ds.coverage.unwrap();
    // let lens: Vec<_> = ds.raw_reads.iter().map(|read| read.seq().len()).collect();
    // let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    // graph.remove_lightweight_edges(1);
    // // First clean up.
    // graph.remove_zero_copy_elements(coverage, &lens, 0.51);
    // graph.z_edge_selection();
    // graph.remove_tips(0.5, 5);
    // graph.zip_up_overclustering();
    // graph.resolve_repeats(&reads, &config, 8f64);
    // // Second clean up.
    // graph.assign_copy_number(coverage, &lens);
    // graph.resolve_repeats(&reads, &config, 4f64);
    // // Third clustering.
    // graph.assign_copy_number(coverage, &lens);
    // graph.z_edge_selection();
    // graph.zip_up_overclustering();
    // graph.resolve_repeats(&reads, &config, 4f64);
    // // Assemble
    // let gfa = assemble_as_gfa(&graph, &config);
    // println!("{}", gfa);
    Ok(())
}

// #[allow(dead_code)]
// fn assemble_as_gfa(
//     graph: &haplotyper::assemble::ditch_graph::DitchGraph,
//     c: &haplotyper::AssembleConfig,
// ) -> gfa::GFA {
//     let header = gfa::Content::Header(gfa::Header::default());
//     let header = gfa::Record::from_contents(header, vec![]);
//     let (segments, edge, group, summaries) = graph.spell(c, 0);
//     for summary in summaries.iter() {
//         let ids: Vec<_> = summary
//             .summary
//             .iter()
//             .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
//             .collect();
//         log::debug!("{}\t{}", summary.id, ids.join("\t"));
//     }
//     let nodes = segments.into_iter().map(|node| {
//         let tags = summaries
//             .iter()
//             .find(|x| x.id == node.sid)
//             .map(|contigsummary| {
//                 let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
//                 let coverage =
//                     gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
//                 let (cp, cpnum) = contigsummary
//                     .summary
//                     .iter()
//                     .filter_map(|elm| elm.copy_number)
//                     .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
//                 let mut tags = vec![coverage];
//                 if cpnum != 0 {
//                     tags.push(gfa::SamTag::new(format!("cp:i:{}", cp / cpnum)));
//                 }
//                 tags
//             })
//             .unwrap_or(vec![]);
//         gfa::Record::from_contents(gfa::Content::Seg(node), tags)
//     });
//     let edges = edge
//         .into_iter()
//         .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
//     let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
//     let records: Vec<_> = std::iter::once(group).chain(nodes).chain(edges).collect();
//     let mut header = vec![header];
//     header.extend(records);
//     gfa::GFA::from_records(header)
// }

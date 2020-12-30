pub mod ditch_graph;
pub mod pileup;
pub mod string_graph;
use definitions::*;
// mod naive_consensus;
use bio_utils::lasttab::LastTAB;
use ditch_graph::*;
use gfa::GFA;
use rayon::prelude::*;
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

use std::collections::{HashMap, HashSet};
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
        // <============ String graph =================>
        // let config = &string_graph::DEFAULT_CONFIG;
        // let mut graph = string_graph::StringGraph::from_dataset(&self, config);
        // assert!(graph.sanity_check());
        // debug!("{:?}", graph);
        // graph.transitive_edge_reduction();
        // assert!(graph.sanity_check());
        // debug!("{:?}", graph);
        // graph.simple_path_reduction();
        // assert!(graph.sanity_check());
        // debug!("{:?}", graph);
        // graph.assemble_as_gfa()
        // <============ Ditch graph ======================>
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
                let (nodes, edges, group, _summaries) = assemble(self, cl, c);
                let mut records = vec![];
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
                let (_, edges, _, summaries) = assemble(&self, cl, c);
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
    Vec<gfa::Segment>,
    Vec<gfa::Edge>,
    gfa::Group,
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
    let mut graph = DitchGraph::new(&reads, c);
    graph.resolve_repeats();
    graph.remove_redundant_edges(2);
    graph.remove_tips();
    graph.collapse_buddle(c);
    graph.remove_small_component(5);
    debug!("{}", graph);
    let (segments, edge, group, summaries) = graph.spell(c, cl);
    let segments = if c.to_polish {
        segments
            .iter()
            .map(|segment| {
                let summary = summaries.iter().find(|s| s.id == segment.sid).unwrap();
                let segment = polish_segment(ds, segment, summary, c);
                let segment = polish_segment(ds, &segment, summary, c);
                correct_short_indels(ds, &segment, summary, c)
            })
            .collect()
    } else {
        segments
    };
    (segments, edge, group, summaries)
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

fn correct_short_indels(
    ds: &DataSet,
    segment: &gfa::Segment,
    summary: &ContigSummary,
    c: &AssembleConfig,
) -> gfa::Segment {
    let reads = get_reads_in_cluster(ds, summary);
    let alignments = match align_reads(segment, &reads, c) {
        Ok(res) => res,
        Err(why) => panic!("{:?}", why),
    };
    use pileup::*;
    let seq = Pileups::convert_into_pileup(&alignments, segment, &reads, c).generate();
    let seq = String::from_utf8(seq).unwrap();
    gfa::Segment::from(segment.sid.clone(), seq.len(), Some(seq))
}

fn align_reads(
    segment: &gfa::Segment,
    reads: &[&RawRead],
    c: &AssembleConfig,
) -> std::io::Result<Vec<LastTAB>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut c_dir = std::env::current_dir()?;
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir)?;
    // Create reference and reads.
    let (reference, reads) = {
        // use bio_utils::fasta;
        let mut reference = c_dir.clone();
        reference.push("segment.fa");
        use std::io::{BufWriter, Write};
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        let mut seq = segment.sequence.as_ref().unwrap().as_bytes().to_vec();
        // TODO: Should be tuned.
        let repeat_masking_config = crate::repeat_masking::RepeatMaskConfig::new(15, 0.0002, 20);
        crate::repeat_masking::mask_repeat_in_seq(&mut seq, &repeat_masking_config);
        let seq = String::from_utf8_lossy(&seq);
        writeln!(&mut wtr, ">{}\n{}", &segment.sid, seq)?;
        // let record = fasta::Record::with_data(&segment.sid, &None, seq);
        // wtr.write_record(&record)?;
        let mut reads_dir = c_dir.clone();
        reads_dir.push("reads.fa");
        // let mut wtr = fasta::Writer::new(std::fs::File::create(&reads_dir)?);
        let mut wtr = std::fs::File::create(&reads_dir).map(BufWriter::new)?;
        let mut read_seqs: Vec<Vec<u8>> = reads.iter().map(|r| r.seq.as_bytes().to_vec()).collect();
        // TODO: Should be tuned.
        let repeat_masking_config = crate::repeat_masking::RepeatMaskConfig::new(15, 0.0002, 200);
        crate::repeat_masking::mask_repeats_in_reads(&mut read_seqs, &repeat_masking_config);
        for (read, seq) in reads.iter().zip(read_seqs) {
            let seq = String::from_utf8_lossy(&seq);
            writeln!(&mut wtr, ">{}\n{}", read.name, seq)?;
            // let id = read.name.to_string();
            // let record = fasta::Record::with_data(&id, &None, read.seq.as_bytes());
            // wtr.write_record(&record)?;
        }
        let reference = reference.into_os_string().into_string().unwrap();
        let reads = reads_dir.into_os_string().into_string().unwrap();
        (reference, reads)
    };
    let alignment = invoke_last(&reference, &reads, &c_dir, c.threads);
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    debug!("Alignment done");
    alignment
}

pub fn invoke_last(
    reference: &str,
    reads: &str,
    c_dir: &std::path::PathBuf,
    threads: usize,
) -> std::io::Result<Vec<LastTAB>> {
    use std::io::{BufWriter, Write};
    let mut db_name = c_dir.clone();
    db_name.push("reference");
    let db_name = db_name.into_os_string().into_string().unwrap();
    // Create database - train - align
    let lastdb = std::process::Command::new("lastdb")
        .args(&["-R", "11", "-Q", "0", &db_name, &reference])
        .output()?;
    if !lastdb.status.success() {
        panic!("lastdb-{}", String::from_utf8_lossy(&lastdb.stderr));
    }
    let p = format!("{}", threads);
    let last_train = std::process::Command::new("last-train")
        .args(&["-P", &p, "-Q", "0", &db_name, &reads])
        .output()
        .unwrap();
    if !last_train.status.success() {
        panic!("last-train-{}", String::from_utf8_lossy(&last_train.stderr));
    }
    let param = {
        let mut param = c_dir.clone();
        param.push("param.par");
        let mut wtr = BufWriter::new(std::fs::File::create(&param).unwrap());
        wtr.write_all(&last_train.stdout).unwrap();
        wtr.flush().unwrap();
        param.into_os_string().into_string().unwrap()
    };
    let lastal = std::process::Command::new("lastal")
        .args(&[
            "-P", &p, "-R", "11", "-Q", "0", "-p", &param, &db_name, &reads,
        ])
        .stdout(std::process::Stdio::piped())
        .spawn();
    let lastal = match lastal {
        Ok(res) => match res.stdout {
            Some(res) => res,
            None => panic!("lastal invoke, but stdout is None."),
        },
        Err(why) => panic!("lastal{:?}", why),
    };
    let last_split = std::process::Command::new("last-split")
        .stdin(std::process::Stdio::from(lastal))
        .stdout(std::process::Stdio::piped())
        .spawn();
    let last_split = match last_split {
        Ok(res) => match res.stdout {
            Some(res) => res,
            None => panic!("last_split is invoked, but the stdout is None"),
        },
        Err(why) => panic!("last_split {:?}", why),
    };
    let maf_convert = match std::process::Command::new("maf-convert")
        .args(&["tab", "-j", "1000"])
        .stdin(std::process::Stdio::from(last_split))
        .output()
    {
        Ok(res) if res.status.success() => res,
        Ok(res) => panic!("maf-convert exit with status {}", res.status),
        Err(why) => panic!("maf_convert {:?}", why),
    };
    let alignments: Vec<_> = String::from_utf8_lossy(&maf_convert.stdout)
        .lines()
        .filter(|e| !e.starts_with('#'))
        .filter_map(|e| LastTAB::from_line(&e))
        .collect();
    Ok(alignments)
}

pub fn polish_by_chunking(
    alignments: &[LastTAB],
    segment: &gfa::Segment,
    reads: &[&RawRead],
    c: &AssembleConfig,
) -> Vec<u8> {
    let len = c.window_size;
    debug!("Recording {} alignments...", alignments.len());
    let segment = segment.sequence.as_ref().unwrap().as_bytes();
    let alignments: HashMap<_, Vec<_>> = {
        let mut result: HashMap<_, Vec<_>> = HashMap::new();
        for aln in alignments
            .into_iter()
            .filter(|aln| aln.seq1_matchlen() > crate::encode::MARGIN)
        {
            result
                .entry(aln.seq2_name().to_string())
                .or_default()
                .push(aln);
        }
        result
    };
    let mut chunks: Vec<Vec<_>> = vec![vec![]; segment.len() / len + 1];
    for read in reads.iter() {
        if let Some(alns) = alignments.get(&read.name) {
            let cs = into_chunks(alns, read, segment, len);
            for (pos, chunk) in cs {
                chunks[pos].push(chunk);
            }
        }
    }
    chunks
        .into_par_iter()
        .enumerate()
        .flat_map(|(idx, mut cs)| {
            if cs.is_empty() {
                return vec![];
            }
            let mut lens: Vec<_> = cs.iter().map(|seq| seq.len()).collect();
            lens.sort();
            let median = lens[lens.len() / 2];
            lens.sort_by_key(|x| x.max(&median) - x.min(&median));
            let mad = median.max(lens[lens.len() / 2]) - median.min(lens[lens.len() / 2]);
            let mad = mad.max(1);
            let min = median.max(mad * 5) - mad * 5;
            let max = median + mad * 5;
            cs.retain(|seq| min < seq.len() && seq.len() < max);
            let template = &segment[idx * len..((idx + 1) * len).min(segment.len())];
            consensus(template, &cs, 10)
        })
        .collect()
}

fn consensus(_template: &[u8], xs: &[Vec<u8>], num: usize) -> Vec<u8> {
    if xs.is_empty() {
        return vec![];
    } else if xs.iter().map(|xs| xs.len()).max().unwrap() < 10 {
        return xs.iter().max_by_key(|x| x.len()).unwrap().to_vec();
    }
    let param = (-2, -2, &|x, y| if x == y { 2 } else { -4 });
    use rand_xoshiro::Xoroshiro128StarStar;
    let cs: Vec<_> = (0..num as u64)
        .into_par_iter()
        .map(|s| {
            use rand::seq::SliceRandom;
            let mut rng: Xoroshiro128StarStar = rand::SeedableRng::seed_from_u64(s);
            let mut cs: Vec<_> = xs.to_vec();
            cs.shuffle(&mut rng);
            let max_len = cs.iter().map(|s| s.len()).max().unwrap_or(0);
            let node_num_thr = (max_len as f64 * 1.5).floor() as usize;
            cs.iter()
                .filter(|c| !c.is_empty())
                //.fold(poa_hmm::POA::new(template, 3.), |x, y| {
                .fold(poa_hmm::POA::default(), |x, y| {
                    let res = if x.nodes().len() > node_num_thr {
                        x.add(y, 1., param).remove_node(0.4)
                    } else {
                        x.add(y, 1., param)
                    };
                    res
                })
                .remove_node(0.3)
                .finalize()
                .consensus()
        })
        .collect();
    let cs: Vec<_> = cs.iter().map(|cs| cs.as_slice()).collect();
    poa_hmm::POA::default()
        .update_thr(&cs, &vec![1.; cs.len()], param, 0.8, 1.5)
        .consensus()
}

// Maybe we can use encoding module to this functionality?
// Or, maybe we just pick the most possible alignment...?
fn into_chunks<'a>(
    alignments: &[&LastTAB],
    read: &'a RawRead,
    segment: &[u8],
    len: usize,
) -> Vec<(usize, Vec<u8>)> {
    // Determine the direction,
    // alignments
    //     .iter()
    //     .flat_map(|aln| {
    //         let read: Vec<_> = if aln.seq2_direction().is_forward() {
    //             read.seq().iter().map(|x| x.to_ascii_uppercase()).collect()
    //         } else {
    //             let seq: Vec<_> = read.seq().iter().map(|x| x.to_ascii_uppercase()).collect();
    //             bio_utils::revcmp(&seq)
    //         };
    //         let query_start = aln.seq2_start();
    //         let ref_start = aln.seq1_start();
    //         use bio_utils::lasttab;
    //         let ops: Vec<_> = aln
    //             .alignment()
    //             .into_iter()
    //             .map(|op| match op {
    //                 lasttab::Op::Seq1In(l) => Op::Ins(*l),
    //                 lasttab::Op::Seq2In(l) => Op::Del(*l),
    //                 lasttab::Op::Match(l) => Op::Match(*l),
    //             })
    //             .collect();
    //         split_reads(&read, query_start, ref_start, ops, len)
    //     })
    //     .collect()
    let read: Vec<_> = read.seq().iter().map(|x| x.to_ascii_uppercase()).collect();
    use crate::encode::join_alignments;
    let merged = [true, false]
        .iter()
        .filter_map(|&direction| {
            let mut alns: Vec<&LastTAB> = alignments
                .iter()
                .copied()
                .filter(|aln| aln.seq2_direction().is_forward() == direction)
                .collect();
            if alns.is_empty() {
                return None;
            }
            let read: Vec<_> = if direction {
                read.clone()
            } else {
                bio_utils::revcmp(&read)
            };
            let (qpos, rpos, ops) = join_alignments(&mut alns, segment, &read);
            let score = ops
                .iter()
                .map(|&e| match e {
                    Op::Del(l) => l as i32 * -2,
                    Op::Ins(l) => l as i32 * -2,
                    Op::Match(l) => l as i32,
                })
                .sum::<i32>();
            Some((qpos, rpos, ops, read, score, direction))
        })
        .max_by_key(|x| x.4);
    if let Some((query_start, refr_start, ops, read, _, _)) = merged {
        split_reads(&read, query_start, refr_start, ops, len)
    } else {
        Vec::with_capacity(0)
    }
}

pub fn split_reads(
    read: &[u8],
    query_start: usize,
    ref_start: usize,
    mut ops: Vec<Op>,
    len: usize,
) -> Vec<(usize, Vec<u8>)> {
    let (mut q_pos, mut r_pos) = (query_start, ref_start);
    let mut target = (ref_start / len + 1) * len;
    ops.reverse();
    let mut chunk_position = vec![query_start];
    let (query_total, refr_total) = ops.iter().fold((q_pos, r_pos), |(q, r), &x| match x {
        Op::Match(l) => (q + l, r + l),
        Op::Del(l) => (q, r + l),
        Op::Ins(l) => (q + l, r),
    });
    loop {
        // Pop until target.
        let last_op = {
            let rest = refr_total - r_pos;
            loop {
                match ops.pop() {
                    Some(Op::Del(l)) => {
                        r_pos += l;
                        if target <= r_pos {
                            break Op::Del(l);
                        }
                    }
                    Some(Op::Ins(l)) => {
                        q_pos += l;
                    }
                    Some(Op::Match(l)) => {
                        r_pos += l;
                        q_pos += l;
                        if target <= r_pos {
                            break Op::Match(l);
                        }
                    }
                    None => {
                        debug!(
                            "ERR\t{}\t{}\t{}\t{}\t{}\t{}",
                            query_total, refr_total, q_pos, r_pos, target, rest
                        );
                        unreachable!()
                    }
                }
            }
        };
        // Push back
        if target < r_pos {
            let overflow = r_pos - target;
            match last_op {
                Op::Del(_) => {
                    ops.push(Op::Del(overflow));
                    r_pos -= overflow;
                }
                Op::Match(_) => {
                    ops.push(Op::Match(overflow));
                    r_pos -= overflow;
                    q_pos -= overflow;
                }
                _ => unreachable!(),
            }
        }
        chunk_position.push(q_pos);
        // Move to next iteration
        let rest = refr_total - r_pos;
        let query_rest = query_total - q_pos;
        if rest < len && query_rest > 0 {
            chunk_position.push(query_total);
            break;
        } else if rest < len && query_rest == 0 {
            break;
        } else {
            target += len;
        }
    }
    let offset = ref_start / len;
    chunk_position
        .windows(2)
        .enumerate()
        .map(|(idx, w)| (idx + offset, read[w[0]..w[1]].to_vec()))
        .collect()
}

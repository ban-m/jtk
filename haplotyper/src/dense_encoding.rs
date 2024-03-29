use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::{assemble::ditch_graph::ContigSummary, ALN_PARAMETER};
const CONS_MIN_LENGTH: usize = 400;
const CONS_MAX_LENGTH: usize = 10_000;
// Directed edge between nodes.
type DEdge = ((u64, u64, bool), (u64, u64, bool));
// (chunk,cluster,direction, if it is `from` part)
type DTip = (u64, u64, bool, bool);
#[derive(Debug, Clone)]
pub struct DenseEncodingConfig {
    pub len: usize,
    pub file: Option<String>,
}
impl DenseEncodingConfig {
    pub fn new(len: usize, file: Option<&str>) -> Self {
        Self {
            len,
            file: file.map(|f| f.to_string()),
        }
    }
}
pub trait DenseEncoding {
    fn dense_encoding(&mut self, config: &DenseEncodingConfig);
}

fn check_raw_read_consistency(ds: &DataSet) {
    let raw_seq: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r)).collect();
    let is_ok = ds.encoded_reads.par_iter().all(|read| {
        let orig = read.recover_raw_read();
        assert!(orig.iter().all(u8::is_ascii_uppercase));
        let raw: Vec<_> = raw_seq[&read.id]
            .seq()
            .iter()
            .map(|x| x.to_ascii_uppercase())
            .collect();
        orig == raw
    });
    assert!(is_ok);
}

const SQUISH_ARI_THR: f64 = 0.5;
const SQUISH_COUNT_THR: usize = 7;
const MAT_ARI: f64 = 4f64;
const MIS_ARI: f64 = -1f64;
impl DenseEncoding for DataSet {
    fn dense_encoding(&mut self, config: &DenseEncodingConfig) {
        check_raw_read_consistency(self);
        self.sanity_check();
        let original_cluster_num: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        let original_assignments = log_original_assignments(self);
        use crate::phmm_likelihood_correction::*;
        use crate::{SquishConfig, SquishErroneousClusters};
        let cor_config = CorrectionConfig::default();
        let squish_config = SquishConfig::new(SQUISH_ARI_THR, SQUISH_COUNT_THR, MAT_ARI, MIS_ARI);
        self.squish_erroneous_clusters(&squish_config);
        self.correct_clustering(&cor_config);
        let new_chunks = encode_polyploid_edges(self, config);
        for read in self.encoded_reads.iter_mut() {
            let orig = &original_assignments[&read.id];
            recover_original_assignments(read, orig);
        }
        self.selected_chunks
            .iter_mut()
            .filter(|c| original_cluster_num.contains_key(&c.id))
            .for_each(|c| c.cluster_num = original_cluster_num[&c.id]);
        let raw_reads: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, r)).collect();
        self.encoded_reads.par_iter_mut().for_each(|read| {
            let mut nodes = Vec::with_capacity(read.nodes.len());
            nodes.append(&mut read.nodes);
            use crate::encode::{nodes_to_encoded_read, remove_overlapping_encoding};
            nodes = remove_overlapping_encoding(nodes);
            let seq = raw_reads[&read.id].seq();
            *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
        });
        self.sanity_check();
        use crate::local_clustering::LocalClustering;
        self.local_clustering_selected(&new_chunks);
    }
}

type TipAndUnit<'a> = HashMap<DTip, Vec<&'a [Chunk]>>;
fn encode_polyploid_edges(ds: &mut DataSet, config: &DenseEncodingConfig) -> HashSet<u64> {
    let edge_chunks = enumerate_polyploid_edges(ds, config);
    let tip_chunks = {
        let mut tip_chunks: TipAndUnit = HashMap::new();
        for (&(from, to), chunks) in edge_chunks.iter() {
            let from = (from.0, from.1, from.2, true);
            let to = (to.0, to.1, to.2, false);
            tip_chunks.entry(from).or_default().push(chunks.as_slice());
            tip_chunks.entry(to).or_default().push(chunks.as_slice());
        }
        tip_chunks
    };
    let readtype = ds.read_type;
    ds.encoded_reads
        .par_iter_mut()
        .for_each(|read| edge_encode(read, &edge_chunks, &tip_chunks, readtype, config));
    let chunk_ids: HashSet<_> = edge_chunks
        .values()
        .flat_map(|x| x.iter().map(|x| x.id))
        .collect();
    {
        let lens: HashMap<_, _> = edge_chunks
            .values()
            .flat_map(|x| x.iter().map(|x| (x.id, x.seq().len())))
            .collect();
        let mut counts: HashMap<_, u32> = lens.keys().map(|&id| (id, 0)).collect();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *counts.entry(node.chunk).or_default() += 1;
        }
        let mut chunk_ids: Vec<_> = chunk_ids.iter().collect();
        chunk_ids.sort();
        for id in chunk_ids.iter() {
            let (count, len) = (counts[id], lens[id]);
            debug!("DE\tCount\t{}\t{}\t{}", id, count, len);
        }
    }
    let mut current_chunks: HashSet<_> = ds.selected_chunks.iter().map(|c| c.id).collect();
    for (_, chunks) in edge_chunks {
        for chunk in chunks {
            if !current_chunks.contains(&chunk.id) {
                current_chunks.insert(chunk.id);
                ds.selected_chunks.push(chunk);
            }
        }
    }
    chunk_ids
}

fn edge_encode(
    read: &mut EncodedRead,
    edges: &EdgeAndUnit,
    tips: &TipAndUnit,
    read_type: definitions::ReadType,
    _config: &DenseEncodingConfig,
) {
    let seq = read.recover_raw_read();
    let inserts = fill_edges_by_new_chunks(read, &seq, edges, tips, &read_type);
    for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
        match idx + accum_inserts {
            pos if pos < read.nodes.len() => read.nodes.insert(idx + accum_inserts, node),
            _ => read.nodes.push(node),
        }
    }
    re_encode_read(read, &seq);
}

fn re_encode_read(read: &mut EncodedRead, seq: &[u8]) {
    if !read.nodes.is_empty() {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use crate::encode::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| n.chunk);
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
}

type NodeLog = (u64, u64, Vec<f64>);
fn log_original_assignments(ds: &DataSet) -> HashMap<u64, Vec<NodeLog>> {
    ds.encoded_reads
        .iter()
        .map(|r| {
            let xs: Vec<_> = r
                .nodes
                .iter()
                .map(|u| (u.chunk, u.cluster, u.posterior.clone()))
                .collect();
            (r.id, xs)
        })
        .collect()
}

// Recover the previous clustering. Note that sometimes the node is added
// so that the length of the read is different from the logged one.
// But it does not removed!
fn recover_original_assignments(read: &mut EncodedRead, log: &[NodeLog]) {
    let mut read = read.nodes.iter_mut();
    'outer: for &(chunk, cluster, ref post) in log {
        for node in &mut read {
            if node.chunk == chunk {
                node.cluster = cluster;
                node.posterior.clear();
                node.posterior.extend(post);
                continue 'outer;
            }
        }
        panic!("Can not find the logged chunks.");
    }
}

pub fn fill_edges_by_new_chunks(
    read: &EncodedRead,
    seq: &[u8],
    edges: &EdgeAndUnit,
    tips: &TipAndUnit,
    read_type: &definitions::ReadType,
) -> Vec<(usize, Node)> {
    const MARGIN: usize = 25;
    let len = seq.len();
    let mut inserts = vec![];
    // Head tip.
    if let Some(head) = read.nodes.first() {
        let (start, end) = (0, (head.position_from_start + MARGIN).min(len));
        let key = (head.chunk, head.cluster, head.is_forward, false);
        if let Some(chunks) = tips.get(&key) {
            // --Tip--|Node[0]>|------
            // --Unit-|ToNode|-----
            for chunk_info in chunks {
                let nodes = encode_edge(seq, start, end, true, chunk_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (0, x)));
            }
        }
        // Here is a bug
        let key = (head.chunk, head.cluster, !head.is_forward, true);
        if let Some(chunks) = tips.get(&key) {
            // |<Node[0]|-Tip--
            // |FromNode|-Unit-
            for chunk_info in chunks {
                let nodes = encode_edge(seq, start, end, false, chunk_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (0, x)));
            }
        }
    }
    for (idx, window) in read.nodes.windows(2).enumerate() {
        let start = window[0].position_from_start + window[0].seq().len();
        let start = (start.max(MARGIN) - MARGIN).min(len);
        let end = window[1].position_from_start;
        let end = (end + MARGIN).min(len);
        let forward = get_forward_d_edge_from_window(window);
        let reverse = get_reverse_d_edge_from_window(window);
        let (chunk_info, direction) = if edges.contains_key(&forward) {
            (&edges[&forward], true)
        } else if edges.contains_key(&reverse) {
            (&edges[&reverse], false)
        } else {
            continue;
        };
        if end <= start {
            let start = window[0].position_from_start;
            let end = window[1].position_from_start;
            let seqlen = window[0].seq().len();
            warn!(
                "TooNear\t{}\t{}\t{}\t{}\t{}",
                read.id, start, end, seqlen, len
            );
            let from = (window[0].chunk, window[0].cluster);
            let to = (window[1].chunk, window[1].cluster);
            warn!("Dump\t{}\t{:?}\t{:?}", read.id, from, to);
            continue;
        }
        let encoded = encode_edge(seq, start, end, direction, chunk_info, read_type);
        for node in encoded {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    // Tail tip
    if let Some(tail) = read.nodes.last() {
        let idx = read.nodes.len();
        let (start, end) = ((tail.position_from_start + tail.seq().len()).min(len), len);
        let (start, end) = (start.saturating_sub(MARGIN), (end + MARGIN).min(len));
        let key = (tail.chunk, tail.cluster, tail.is_forward, true);
        if let Some(chunks) = tips.get(&key) {
            // | Last>  |-Tip--
            // |FromNode|-Unit-
            for chunk_info in chunks {
                let nodes = encode_edge(seq, start, end, true, chunk_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (idx + 1, x)));
            }
        }
        let key = (tail.chunk, tail.cluster, !tail.is_forward, false);
        if let Some(chunks) = tips.get(&key) {
            // --Tip-|<Node[0]|
            // -Unit-|ToNode|
            for chunk_info in chunks {
                let nodes = encode_edge(seq, start, end, false, chunk_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (idx + 1, x)));
            }
        }
    }
    inserts
}

fn write_to_file(
    records: &[gfa::Record],
    summaries: &[ContigSummary],
    de_config: &DenseEncodingConfig,
) {
    if let Some(file) = de_config.file.as_ref() {
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![].into());
        let records = std::iter::once(header)
            .chain(records.iter().cloned())
            .collect();
        let gfa = gfa::GFA::from_records(records);
        if let Ok(mut wtr) = std::fs::File::create(file).map(std::io::BufWriter::new) {
            use std::io::Write;
            if let Err(why) = writeln!(wtr, "{}", gfa) {
                eprintln!("{:?}", why);
            }
        }
    }
    for summary in summaries {
        let (copy_num, tig_num) = summary
            .summary
            .iter()
            .filter_map(|s| s.copy_number)
            .fold((0, 0), |(c, x), copynum| (c + copynum, x + 1));
        let copy_num = match tig_num {
            0 => 0,
            _ => (copy_num as f64 / tig_num as f64).round() as usize,
        };
        let ids: Vec<_> = summary
            .summary
            .iter()
            .map(|elm| format!("{}-{}", elm.chunk, elm.cluster))
            .collect();
        debug!("DRAFT2\t{}\t{}\t{}", summary.id, copy_num, ids.join("\t"));
    }
}

type EdgeAndUnit = HashMap<DEdge, Vec<Chunk>>;
fn enumerate_polyploid_edges(ds: &DataSet, de_config: &DenseEncodingConfig) -> EdgeAndUnit {
    use crate::assemble::*;
    let msr = ds.read_type.weak_span_reads();
    let min_lk = ds.read_type.weak_llr();
    let config = AssembleConfig::new(1000, false, true, msr, min_lk, false, None);
    let (records, summaries) = assemble(ds, &config);
    write_to_file(&records, &summaries, de_config);
    let multicopy_contigs: HashMap<_, _> = summaries
        .iter()
        .filter(|summary| !summary.summary.is_empty())
        .filter(|summary| {
            let mut has_edges = [false, false];
            for rec in records.iter() {
                if let gfa::Content::Edge(ref edge) = &rec.content {
                    if edge.sid1.id == summary.id {
                        has_edges[edge.end1.is_last as usize] |= true;
                    }
                    if edge.sid2.id == summary.id {
                        has_edges[edge.end2.is_last as usize] |= true;
                    }
                }
            }
            has_edges[0] && has_edges[1]
        })
        .filter_map(|summary| {
            let (total_cp, num) = summary
                .summary
                .iter()
                .filter_map(|x| x.copy_number)
                .fold((0, 0), |(cp, n), x| (cp + x, n + 1));
            let copy_num = total_cp / num;
            (1 < copy_num).then(|| (summary.id.clone(), copy_num))
        })
        .collect();
    let edges: HashMap<_, _> = summaries
        .iter()
        .flat_map(|summary| match multicopy_contigs.get(&summary.id) {
            Some(&copy_number) => summary
                .summary
                .windows(2)
                .map(|w| {
                    let from = (w[0].chunk, w[0].cluster, w[0].strand);
                    let to = (w[1].chunk, w[1].cluster, w[1].strand);
                    ((from, to), copy_number)
                })
                .collect(),
            None => Vec::new(),
        })
        .collect();
    debug!("DE\t{}\tEDGES", edges.len());
    let cov_thr = ds.coverage.unwrap().ceil() as usize / 5;
    let mean_chunk_len = {
        let sum: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
        sum / ds.selected_chunks.len()
    };
    let mut newly_defined_chunk = HashMap::new();
    let mut max_chunk_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    for (key, consensus, copy_num) in take_consensus_on_multitig(ds, &edges, cov_thr) {
        let chunk_num = (consensus.len() as f64 / mean_chunk_len as f64).ceil();
        let chunk_len = (consensus.len() as f64 / chunk_num).ceil() as usize;
        let chunks: Vec<_> = consensus
            .chunks(chunk_len)
            .map(|seq| {
                max_chunk_id += 1;
                Chunk::new(max_chunk_id, seq.to_vec(), copy_num)
            })
            .collect();
        let edge = format!("({},{})-({},{})", key.0 .0, key.0 .2, key.1 .0, key.1 .2);
        let len = consensus.len();
        let (start, end) = (max_chunk_id - chunks.len() as u64, max_chunk_id);
        debug!("DE\tInNode\t{len}\t{start}\t{end}\t{edge}",);
        newly_defined_chunk.insert(key, chunks);
    }
    let consensi =
        take_consensus_to_multitig(ds, &records, &summaries, &multicopy_contigs, cov_thr);
    for (edges, consensus, copy_num) in consensi {
        let chunk_num = (consensus.len() as f64 / mean_chunk_len as f64).ceil();
        let chunk_len = (consensus.len() as f64 / chunk_num).ceil() as usize;
        let chunks: Vec<_> = consensus
            .chunks(chunk_len)
            .map(|seq| {
                max_chunk_id += 1;
                Chunk::new(max_chunk_id, seq.to_vec(), copy_num)
            })
            .collect();
        let len = consensus.len();
        let (start, end) = (max_chunk_id - chunks.len() as u64, max_chunk_id);
        for key in edges.iter() {
            let &((u1, c1, d1), (u2, c2, d2)) = key;
            let edge = format!("({u1},{c1},{d1})-({u2},{c2},{d2})");
            debug!("DE\tBetEdge\t{len}\t{start}\t{end}\t{edge}",);
            newly_defined_chunk.insert(*key, chunks.clone());
        }
    }
    debug!("DE\tDefined\t{}", newly_defined_chunk.len());
    newly_defined_chunk
}

fn take_consensus_on_multitig(
    ds: &DataSet,
    edges: &HashMap<DEdge, usize>,
    cov_thr: usize,
) -> Vec<(DEdge, Vec<u8>, usize)> {
    let mut consensi_materials: HashMap<_, Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        assert_eq!(read.nodes.len(), read.edges.len() + 1);
        for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
            assert_eq!((edge.from, edge.to), (w[0].chunk, w[1].chunk));
            let forward = get_forward_d_edge_from_window(w);
            let reverse = get_reverse_d_edge_from_window(w);
            if edges.contains_key(&forward) {
                let label = edge.label().to_vec();
                consensi_materials.entry(forward).or_default().push(label);
            } else if edges.contains_key(&reverse) {
                let label = bio_utils::revcmp(edge.label());
                consensi_materials.entry(reverse).or_default().push(label);
            }
        }
    }
    debug!("DE\tCand\t{}", consensi_materials.len());
    let mut consensi: Vec<_> = consensi_materials
        .into_par_iter()
        .filter_map(|(key, mut seqs)| {
            seqs.iter_mut()
                .for_each(|xs| xs.iter_mut().for_each(u8::make_ascii_uppercase));
            let mean = seqs.iter().map(|x| x.len()).sum::<usize>() / seqs.len();
            let len = seqs.len();
            let cp = *edges.get(&key).unwrap();
            let cons = consensus(seqs, cov_thr).map(|s| (key, s, cp));
            if cons.is_none() {
                debug!("DE\tEDGE\t{}\t{}\tNG", len, mean);
            } else {
                debug!("DE\tEDGE\t{}\t{}\tOK", len, mean);
            }
            cons
        })
        .collect();
    debug!("DE\tCons\t{}", consensi.len());
    consensi.sort_unstable_by_key(|x| x.0);
    consensi
}

fn take_consensus_to_multitig(
    ds: &DataSet,
    records: &[gfa::Record],
    summaries: &[ContigSummary],
    multicopy_contigs: &HashMap<String, usize>,
    cov_thr: usize,
) -> Vec<(Vec<DEdge>, Vec<u8>, usize)> {
    fn find(
        summaries: &[ContigSummary],
        refid: &String,
        pos: gfa::Position,
        is_from: bool,
    ) -> (u64, u64, bool) {
        let summary = summaries.iter().find(|s| &s.id == refid).unwrap();
        let node = match pos.is_last {
            true => summary.summary.last().unwrap(),
            false => summary.summary.first().unwrap(),
        };
        match (pos.is_last, is_from) {
            (true, true) => (node.chunk, node.cluster, node.strand),
            (false, true) => (node.chunk, node.cluster, !node.strand),
            (true, false) => (node.chunk, node.cluster, !node.strand),
            (false, false) => (node.chunk, node.cluster, node.strand),
        }
    }
    let mut into_multitig_edges: HashMap<_, Vec<_>> = HashMap::new();
    for rec in records.iter() {
        if let gfa::Content::Edge(ref edge) = &rec.content {
            if let Some(cp) = multicopy_contigs.get(&edge.sid1.id) {
                let from = find(summaries, &edge.sid1.id, edge.end1, true);
                let to = find(summaries, &edge.sid2.id, edge.end2, false);
                into_multitig_edges.entry(from).or_default().push((to, *cp))
            }
            if let Some(cp) = multicopy_contigs.get(&edge.sid2.id) {
                let from = find(summaries, &edge.sid2.id, edge.end2, true);
                let to = find(summaries, &edge.sid1.id, edge.end1, false);
                into_multitig_edges.entry(from).or_default().push((to, *cp));
            }
        }
    }
    into_multitig_edges.retain(|_, tos| {
        let to = tos.first().map(|((u, _, s), cp)| (*u, *s, *cp)).unwrap();
        tos.iter().all(|((u, _, s), cp)| (*u, *s, *cp) == to)
    });
    into_multitig_edges
        .into_par_iter()
        .filter_map(|(from, tos)| {
            let cp = tos[0].1;
            let edges: Vec<_> = tos.iter().map(|&(to, _)| (from, to)).collect();
            let mut labels = vec![];
            for read in ds.encoded_reads.iter() {
                assert_eq!(read.nodes.len(), read.edges.len() + 1);
                for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
                    assert_eq!((edge.from, edge.to), (w[0].chunk, w[1].chunk));
                    let forward = get_forward_d_edge_from_window(w);
                    let reverse = get_reverse_d_edge_from_window(w);
                    if edges.iter().any(|&e| e == forward) {
                        let lab: Vec<_> = edge.label().iter().map(u8::to_ascii_uppercase).collect();
                        labels.push(lab);
                    }
                    if edges.iter().any(|&e| e == reverse) {
                        let mut lab = bio_utils::revcmp(edge.label());
                        lab.iter_mut().for_each(u8::make_ascii_uppercase);
                        labels.push(lab);
                    }
                }
            }
            consensus(labels, cov_thr).map(|cons| (edges, cons, cp))
        })
        .collect()
}

fn consensus(mut seqs: Vec<Vec<u8>>, cov_thr: usize) -> Option<Vec<u8>> {
    let pos = seqs.len() / 2;
    let median = seqs.select_nth_unstable_by_key(pos, |x| x.len()).1.len();
    let (upper, lower) = (2 * median, median.max(CONS_MIN_LENGTH) / 2);
    if upper <= lower || CONS_MAX_LENGTH < median {
        return None;
    }
    let idx = seqs.iter().position(|x| x.len() == median).unwrap();
    seqs.swap(0, idx);
    seqs.retain(|x| (lower..upper).contains(&x.len()));
    if seqs.len() <= cov_thr {
        return None;
    }
    let draft = seqs[0].to_vec();
    // let draft = kiley::polish_by_pileup(&seqs[0], &seqs[1..]);
    let task = edlib_sys::AlignTask::Alignment;
    let mode = edlib_sys::AlignMode::Global;
    let mut ops: Vec<_> = seqs
        .iter()
        .map(|x| {
            let aln = edlib_sys::align(x, &draft, mode, task);
            crate::misc::edlib_to_kiley(aln.operations().unwrap())
        })
        .collect();
    let mean_len = seqs.iter().map(|x| x.len()).sum::<usize>() / seqs.len();
    let band_width = (mean_len / 20).max(10).min(50);
    let draft =
        kiley::bialignment::guided::polish_until_converge_with(&draft, &seqs, &mut ops, band_width);
    let cons =
        kiley::bialignment::guided::polish_until_converge_with(&draft, &seqs, &mut ops, band_width);
    (CONS_MIN_LENGTH < cons.len()).then_some(cons)
}

// w: windows of nodes with 2 length.
fn get_forward_d_edge_from_window(w: &[Node]) -> DEdge {
    let from = (w[0].chunk, w[0].cluster, w[0].is_forward);
    let to = (w[1].chunk, w[1].cluster, w[1].is_forward);
    (from, to)
}

// w: windows of nodes with 2 length.
fn get_reverse_d_edge_from_window(w: &[Node]) -> DEdge {
    let from = (w[1].chunk, w[1].cluster, !w[1].is_forward);
    let to = (w[0].chunk, w[0].cluster, !w[0].is_forward);
    (from, to)
}

fn merge_chunks(chunks: &[Chunk]) -> (Vec<u8>, Vec<usize>) {
    let contig: Vec<_> = chunks.iter().map(|x| x.seq()).fold(Vec::new(), |mut x, y| {
        x.extend(y);
        x
    });
    let (break_points, _): (Vec<_>, _) =
        chunks
            .iter()
            .map(|x| x.seq().len())
            .fold((Vec::new(), 0), |(mut acc, x), y| {
                acc.push(x + y);
                (acc, x + y)
            });
    (contig, break_points)
}

fn remove_leading_insertions(ops: &mut Vec<kiley::Op>) -> usize {
    let mut head_ins = 0;
    ops.reverse();
    while let Some(&kiley::Op::Ins) = ops.last() {
        ops.pop().unwrap();
        head_ins += 1;
    }
    ops.reverse();
    head_ins
}

fn alignment_identity(ops: &[kiley::Op]) -> f64 {
    let mat = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
    mat as f64 / ops.len() as f64
}
// Note that the seq[start..end] can be much longer than the contig itself....
fn encode_edge(
    seq: &[u8],
    start: usize,
    end: usize,
    is_forward: bool,
    chunks: &[Chunk],
    read_type: &definitions::ReadType,
) -> Vec<definitions::Node> {
    let (contig, break_points) = merge_chunks(chunks);
    // seq is equal to seq[start..end], revcmped if is_forward is false.
    let band = read_type.band_width(contig.len());
    let ((start, end, seq), (ctg_start, ctg_end, _), mut ops) =
        tune_position(start, end, seq, is_forward, &contig, band);
    let mut xpos = ctg_start;
    // Nearest break;
    let mut target_idx = match break_points.iter().position(|&x| ctg_start < x) {
        Some(i) => i,
        None => return Vec::new(),
    };
    // Current edit operations inside the focal chunk.
    let mut alignments = match target_idx {
        0 => vec![kiley::Op::Del; ctg_start],
        i => vec![kiley::Op::Del; ctg_start - break_points[i - 1]],
    };
    // Push deletion operations up to the last base.
    ops.extend(std::iter::repeat(kiley::Op::Del).take(contig.len() - ctg_end));
    // Split the alignment into encoded nodes.
    // Current position of the query.
    let mut ypos = remove_leading_insertions(&mut ops);
    // Encoded nodes.
    let mut nodes = vec![];
    let sim_thr = read_type.sim_thr();
    for op in ops {
        match op {
            kiley::Op::Match | kiley::Op::Mismatch => {
                xpos += 1;
                ypos += 1;
            }
            kiley::Op::Del => xpos += 1,
            kiley::Op::Ins => ypos += 1,
        }
        alignments.push(op);
        if xpos == break_points[target_idx] {
            // Reached the boundary.
            if 1f64 - sim_thr < alignment_identity(&alignments) {
                let chunk = &chunks[target_idx];
                let ylen = alignments.iter().filter(|&&x| x != kiley::Op::Del).count();
                let position_from_start = match is_forward {
                    true => start + ypos - ylen,
                    false => end - ypos,
                };
                let seq = seq[ypos - ylen..ypos].to_vec();
                let cl = chunk.cluster_num;
                let cigar = crate::misc::kiley_op_to_ops(&alignments).0;
                let node = Node::new(chunk.id, is_forward, seq, cigar, position_from_start, cl);
                nodes.push(node);
            }
            // Refresh.
            target_idx += 1;
            alignments.clear();
        }
        if target_idx == chunks.len() {
            // This is needed, as sometimes only insertions would be remain.
            break;
        }
    }
    if !is_forward {
        nodes.reverse()
    }
    nodes
}

// Note that the seq and the contig is sometimes like:
// seq:----------------------
// ctg:    --------------
// And sometimes like:
// seq:----------(Tip)
// ctg:   -------------------(Long edge)
// This.
// So, we need two-round infix alignment so that, in first alignment, the sequence would be truncted,
// and the second round the contig would be tructed.
type TunedPosition<'a> = (
    (usize, usize, Vec<u8>),
    (usize, usize, &'a [u8]),
    Vec<kiley::Op>,
);
fn tune_position<'a>(
    start: usize,
    end: usize,
    seq: &[u8],
    is_forward: bool,
    contig: &'a [u8],
    band: usize,
) -> TunedPosition<'a> {
    assert!(start < end);
    let mut seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    // First, let's truncate the `seq`
    let alignment = edlib_sys::align(contig, &seq, mode, task);
    let (seq_start, seq_end) = alignment.location().unwrap();
    let seq_end = seq_end + 1;
    // Modefy seq into truncate version.
    let pop_num = seq.len() - seq_end;
    (0..pop_num).map(|_| seq.pop().unwrap()).count();
    seq.reverse();
    (0..seq_start).map(|_| seq.pop().unwrap()).count();
    seq.reverse();
    // In the original coordinate.
    let (seq_start, seq_end) = match is_forward {
        true => (start + seq_start, start + seq_end),
        false => (end - seq_end, end - seq_start),
    };
    // Second, let's truncate the contig.
    let alignment = edlib_sys::align(&seq, contig, mode, task);
    let (ctg_start, ctg_end) = alignment.location().unwrap();
    let ctg_end = ctg_end + 1;
    let contig = &contig[ctg_start..ctg_end];
    let alignment = alignment.operations().unwrap();
    let edlib_to_op = {
        use kiley::Op::*;
        [Match, Ins, Del, Mismatch]
    };
    let ops: Vec<_> = alignment.iter().map(|&x| edlib_to_op[x as usize]).collect();
    let (_, ops) =
        kiley::bialignment::guided::global_guided(contig, &seq, &ops, band, ALN_PARAMETER);
    ((seq_start, seq_end, seq), (ctg_start, ctg_end, contig), ops)
}

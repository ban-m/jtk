use rand::{prelude::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128PlusPlus;

use crate::consensus::Alignment;
use std::{collections::HashMap, io::BufRead};

pub fn polish_segmnents(
    reads: &str,
    contig: &str,
    alignments: &str,
    format: &str,
    window_size: usize,
    seed: u64,
) -> std::io::Result<()> {
    let reads = parse_input(reads)?;
    let reads_id: HashMap<_, _> = reads
        .iter()
        .enumerate()
        .map(|(i, (id, _))| (id.clone(), i as u64))
        .collect();
    let contigs = Contig::parse(contig)?;
    let segments = contigs.map();
    let mut alignments = match format {
        "sam" => parse_sam(&reads, &reads_id, &segments, alignments, seed)?,
        "paf" => parse_paf(&reads, &reads_id, &segments, alignments, seed)?,
        _ => panic!("{} is not a valid format", format),
    };
    let segments: Vec<_> = alignments
        .iter_mut()
        .map(|(id, alns)| {
            let segment = segments.get(id).unwrap();
            let polished = polish_segment(id, segment, alns, window_size, seed);
            (id.clone(), polished)
        })
        .collect();
    contigs.output(segments);
    Ok(())
}

fn polish_segment(
    sid: &str,
    draft: &[u8],
    alns: &mut [Alignment],
    window_size: usize,
    seed: u64,
) -> Vec<u8> {
    let config = crate::consensus::PolishConfig::new(seed, 2, 50, window_size, 50, 2);
    let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let models = kiley::hmm::guided::PairHiddenMarkovModelOnStrands::new(hmm.clone(), hmm);
    warn!("TODO: Train parameters.");
    crate::consensus::polish(sid, draft, alns, &models, &config)
}

type ReadToId = HashMap<String, u64>;
type Read = (String, Vec<u8>);
type AlignmentBucket = HashMap<String, Vec<Alignment>>;
fn parse_sam(
    reads: &[Read],
    readsid: &ReadToId,
    segments: &HashMap<String, &[u8]>,
    alignments: &str,
    seed: u64,
) -> std::io::Result<AlignmentBucket> {
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    let mut alignments: Vec<_> = if alignments == "-" {
        let stdio = std::io::stdin();
        std::io::BufReader::new(stdio.lock())
            .lines()
            .filter_map(|l| l.ok())
            .filter_map(|l| bio_utils::sam::Sam::new(&l))
            .collect()
    } else {
        std::fs::File::open(alignments)
            .map(std::io::BufReader::new)?
            .lines()
            .filter_map(|l| l.ok())
            .filter_map(|l| bio_utils::sam::Sam::new(&l))
            .collect()
    };
    alignments.retain(|r| r.r_name() != "*");
    alignments.sort_by(|r1, r2| r1.q_name().cmp(r2.q_name()));
    let mut buffer: Vec<&bio_utils::sam::Sam> = vec![];
    let mut prev = None;
    let mut alignments_on_contigs: HashMap<_, _> =
        segments.keys().map(|id| (id.clone(), vec![])).collect();
    for aln in alignments.iter() {
        if prev != Some(aln.q_name()) && !buffer.is_empty() {
            let q_name = buffer[0].q_name();
            let id = readsid[q_name];
            let read = &reads[id as usize];
            register_sam(&mut alignments_on_contigs, &mut buffer, id, read, &mut rng);
        }
        buffer.push(aln);
        prev = Some(aln.q_name());
    }
    Ok(alignments_on_contigs)
}

fn register_sam<R: Rng>(
    alignments: &mut AlignmentBucket,
    buffer: &mut Vec<&bio_utils::sam::Sam>,
    id: u64,
    read: &Read,
    rng: &mut R,
) {
    let len = read.1.len();
    while !buffer.is_empty() {
        let choices: Vec<_> = (0..buffer.len()).collect();
        let scores: Vec<_> = buffer
            .iter()
            .map(|aln| -> i64 {
                aln.cigar()
                    .iter()
                    .map(|op| match *op {
                        bio_utils::sam::Op::Insertion(l)
                        | bio_utils::sam::Op::Deletion(l)
                        | bio_utils::sam::Op::Mismatch(l) => -(l as i64),
                        bio_utils::sam::Op::Align(l) | bio_utils::sam::Op::Match(l) => l as i64,
                        _ => 0,
                    })
                    .sum()
            })
            .collect();
        let max = scores.iter().max().unwrap();
        let idx = *choices
            .choose_weighted(rng, |&k| ((scores[k] - max) as f64).exp())
            .unwrap();
        let aln = buffer.remove(idx);
        buffer.retain(|a| overlap_frac_sam(aln, a, len) < 0.1);
        alignments
            .get_mut(aln.r_name())
            .unwrap()
            .push(sam_to_aln(id, read.1.as_slice(), aln));
    }
}

fn overlap_frac_sam(aln1: &bio_utils::sam::Sam, aln2: &bio_utils::sam::Sam, len: usize) -> f64 {
    let (s1, e1) = match (aln1.mapped_region(), aln1.is_forward()) {
        ((start, end), true) => (start, end),
        ((start, end), false) => (len - end, len - start),
    };
    let (s2, e2) = match (aln2.mapped_region(), aln2.is_forward()) {
        ((start, end), true) => (start, end),
        ((start, end), false) => (len - end, len - start),
    };
    let aln_size = e2 - s2;
    ((e1.min(e2).saturating_sub(s1.max(s2))) as f64) / aln_size as f64
}

fn sam_to_aln(id: u64, read: &[u8], aln: &bio_utils::sam::Sam) -> Alignment {
    let contig = aln.r_name().to_string();
    let ctg_range = aln.get_range();
    let (start, end) = aln.mapped_region();
    let len = read.len();
    let query = match aln.is_forward() {
        true => read[start..end].to_vec(),
        false => bio_utils::revcmp(&read[len - end..len - start]),
    };
    let is_forward = aln.is_forward();
    let mut ops = vec![];
    let mut saw_head_clip = false;
    let (mut head_clip, mut tail_clip) = (0, 0);
    for op in aln.cigar() {
        use bio_utils::sam;
        match op {
            sam::Op::Insertion(l) => ops.extend(std::iter::repeat(kiley::Op::Ins).take(l)),
            sam::Op::Deletion(l) => ops.extend(std::iter::repeat(kiley::Op::Del).take(l)),
            sam::Op::Align(l) | sam::Op::Match(l) | sam::Op::Mismatch(l) => {
                ops.extend(std::iter::repeat(kiley::Op::Match).take(l))
            }
            sam::Op::HardClip(l) | sam::Op::SoftClip(l) if !saw_head_clip => {
                saw_head_clip = true;
                head_clip = l;
            }
            sam::Op::HardClip(l) | sam::Op::SoftClip(l) => tail_clip = l,
            _ => {}
        }
    }
    let clips = (head_clip, tail_clip);
    Alignment::new(id, contig, ctg_range, clips, query, ops, is_forward)
}

fn parse_paf(
    reads: &[Read],
    readsid: &ReadToId,
    segments: &HashMap<String, &[u8]>,
    alignments: &str,
    seed: u64,
) -> std::io::Result<AlignmentBucket> {
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    let mut alignments: Vec<_> = if alignments == "-" {
        let stdio = std::io::stdin();
        std::io::BufReader::new(stdio.lock())
            .lines()
            .filter_map(|l| l.ok())
            .filter_map(|l| bio_utils::paf::PAF::new(&l))
            .collect()
    } else {
        std::fs::File::open(alignments)
            .map(std::io::BufReader::new)?
            .lines()
            .filter_map(|l| l.ok())
            .filter_map(|l| bio_utils::paf::PAF::new(&l))
            .collect()
    };
    alignments.retain(|r| r.tname != "*");
    alignments.sort_by(|r1, r2| r1.qname.cmp(&r2.qname));
    let mut buffer: Vec<&bio_utils::paf::PAF> = vec![];
    let mut prev = None;
    let mut alignments_on_contigs: HashMap<_, _> =
        segments.keys().map(|id| (id.clone(), vec![])).collect();
    for aln in alignments.iter() {
        if prev != Some(&aln.qname) && !buffer.is_empty() {
            let q_name = &buffer[0].qname;
            let id = readsid[q_name];
            let read = &reads[id as usize];
            register_paf(&mut alignments_on_contigs, &mut buffer, id, read, &mut rng);
        }
        buffer.push(aln);
        prev = Some(&aln.qname);
    }
    Ok(alignments_on_contigs)
}

fn register_paf<R: Rng>(
    alignments: &mut AlignmentBucket,
    buffer: &mut Vec<&bio_utils::paf::PAF>,
    id: u64,
    read: &Read,
    rng: &mut R,
) {
    while !buffer.is_empty() {
        let choices: Vec<_> = (0..buffer.len()).collect();
        let scores: Vec<_> = buffer
            .iter()
            .map(|aln| -> i64 {
                let cigar = aln.get_tag("cg").unwrap().1;
                bio_utils::sam::parse_cigar_string(cigar)
                    .iter()
                    .map(|op| match *op {
                        bio_utils::sam::Op::Match(l) | bio_utils::sam::Op::Align(l) => l as i64,
                        bio_utils::sam::Op::Insertion(l)
                        | bio_utils::sam::Op::Deletion(l)
                        | bio_utils::sam::Op::Mismatch(l) => -(l as i64),
                        _ => 0,
                    })
                    .sum()
            })
            .collect();
        let max = scores.iter().max().unwrap();
        let idx = *choices
            .choose_weighted(rng, |&k| ((scores[k] - max) as f64).exp())
            .unwrap();
        let aln = buffer.remove(idx);
        buffer.retain(|p| overlap_frac_paf(aln, p) < 0.1);
        alignments
            .get_mut(&aln.tname)
            .unwrap()
            .push(paf_to_aln(id, read.1.as_slice(), aln));
    }
}

fn overlap_frac_paf(aln1: &bio_utils::paf::PAF, aln2: &bio_utils::paf::PAF) -> f64 {
    let (s1, e1) = (aln1.qstart, aln1.qend);
    let (s2, e2) = (aln2.qstart, aln2.qend);
    (e1.min(e2).saturating_sub(s1.max(s2)) as f64) / (e2 - s2) as f64
}

fn paf_to_aln(id: u64, read: &[u8], aln: &bio_utils::paf::PAF) -> Alignment {
    let contig = aln.tname.clone();
    let ctg_range = (aln.tstart, aln.tend);
    let (start, end) = (aln.qstart, aln.qend);
    let query = match aln.relstrand {
        true => read[start..end].to_vec(),
        false => bio_utils::revcmp(&read[start..end]),
    };
    let clips = match aln.relstrand {
        true => (start, read.len() - end),
        false => (read.len() - end, start),
    };
    let is_forward = aln.relstrand;
    let mut ops = vec![];
    let cigar = aln.get_tag("cg").unwrap().1;
    for op in bio_utils::sam::parse_cigar_string(cigar) {
        use bio_utils::sam;
        match op {
            sam::Op::Insertion(l) => ops.extend(std::iter::repeat(kiley::Op::Ins).take(l)),
            sam::Op::Deletion(l) => ops.extend(std::iter::repeat(kiley::Op::Del).take(l)),
            sam::Op::Align(l) | sam::Op::Match(l) | sam::Op::Mismatch(l) => {
                ops.extend(std::iter::repeat(kiley::Op::Match).take(l))
            }
            _ => {}
        }
    }
    Alignment::new(id, contig, ctg_range, clips, query, ops, is_forward)
}

#[derive(Debug, Clone)]
enum Contig {
    Fasta(Vec<(String, Vec<u8>)>),
    Gfa(gfa::GFA),
}

impl Contig {
    fn parse(contig: &str) -> std::io::Result<Self> {
        if contig.ends_with("gfa") {
            let rdr = std::fs::File::open(contig)?;
            Ok(Contig::Gfa(gfa::GFA::from_reader(rdr)))
        } else if contig.ends_with('a') {
            let records: Vec<_> = bio_utils::fasta::parse_into_vec(contig)?
                .into_iter()
                .map(|read| {
                    let (id, _, seq) = read.into();
                    (id, seq.into_bytes())
                })
                .collect();
            Ok(Contig::Fasta(records))
        } else {
            panic!("{} is not a valid extension for contigs", contig);
        }
    }
    fn map(&self) -> HashMap<String, &[u8]> {
        match self {
            Contig::Fasta(records) => records
                .iter()
                .map(|(id, seq)| (id.clone(), seq.as_slice()))
                .collect(),
            Contig::Gfa(gfa) => gfa
                .iter()
                .filter_map(|record| match &record.content {
                    gfa::Content::Seg(seg) => Some(seg),
                    _ => None,
                })
                .filter_map(|seg| {
                    seg.sequence
                        .as_ref()
                        .map(|seq| (seg.sid.clone(), seq.as_bytes()))
                })
                .collect(),
        }
    }
    fn output(self, mut segments: Vec<(String, Vec<u8>)>) {
        match self {
            Contig::Fasta(records) => {
                for (id, _) in records.iter() {
                    let idx = segments.iter().position(|x| id == &x.0).unwrap();
                    let seg = segments.remove(idx);
                    let seg = std::str::from_utf8(&seg.1).unwrap();
                    println!(">{id}\n{seg}");
                }
            }
            Contig::Gfa(mut gfa) => {
                let seg_lens: HashMap<_, _> = segments
                    .iter()
                    .map(|(id, seq)| (id.clone(), seq.len()))
                    .collect();
                for record in gfa.iter_mut() {
                    match &mut record.content {
                        gfa::Content::Seg(seg) => {
                            seg.slen = seg_lens[&seg.sid] as u64;
                            let idx = segments.iter().position(|x| seg.sid == x.0).unwrap();
                            let (_, polished) = segments.remove(idx);
                            seg.sequence = Some(String::from_utf8(polished).unwrap());
                        }
                        gfa::Content::Edge(edge) => {
                            let s1len = seg_lens[&edge.sid1.id];
                            let s2len = seg_lens[&edge.sid2.id];
                            if edge.beg1.is_last || edge.end1.is_last {
                                edge.beg1.pos = s1len;
                                edge.end1.pos = s1len;
                            }
                            if edge.beg2.is_last || edge.end2.is_last {
                                edge.beg2.pos = s2len;
                                edge.end2.pos = s2len;
                            }
                        }
                        _ => {}
                    }
                }
                println!("{gfa}");
            }
        }
    }
}

fn parse_input(reads: &str) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    if reads.ends_with('a') {
        bio_utils::fasta::parse_into_vec(reads).map(|records| {
            records
                .into_iter()
                .map(|read| {
                    let (id, _, seq) = read.into();
                    (id, seq.into_bytes())
                })
                .collect()
        })
    } else if reads.ends_with('q') {
        bio_utils::fastq::parse_into_vec(reads).map(|records| {
            records
                .into_iter()
                .map(|read| {
                    let (id, seq, _) = read.into();
                    (id, seq)
                })
                .collect()
        })
    } else {
        panic!("{} is not a valid read file", reads);
    }
}

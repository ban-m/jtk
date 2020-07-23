use bio_utils::lasttab;
use bio_utils::lasttab::LastTAB;
use std::collections::HashMap;
pub trait Encode {
    fn encode(self, threads: usize) -> Self;
}

impl Encode for definitions::DataSet {
    fn encode(mut self, threads: usize) -> Self {
        let alignments = match last_alignment(&self, threads) {
            Ok(res) => res,
            Err(why) => panic!("{:?}:Encoding step", why),
        };
        let alignments_each_reads: HashMap<String, Vec<&LastTAB>> = distribute(&alignments);
        let encoded_reads: Vec<_> = self
            .raw_reads
            .iter()
            .filter_map(|read| {
                let alns = alignments_each_reads.get(&read.name)?;
                encode(read, alns, &self.selected_chunks)
            })
            .collect();
        debug!("Encoding {} reads.", encoded_reads.len());
        self.encoded_reads = encoded_reads;
        self
    }
}
fn last_alignment(ds: &definitions::DataSet, p: usize) -> std::io::Result<Vec<LastTAB>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 10_000;
    let mut c_dir = std::env::current_dir()?;
    use std::io::{BufWriter, Write};
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir)?;
    // Create reference and reads.
    let (reference, reads) = {
        use bio_utils::fasta;
        let mut reference = c_dir.clone();
        reference.push("units.fa");
        let mut wtr = fasta::Writer::new(std::fs::File::create(&reference)?);
        for unit in ds.selected_chunks.iter() {
            let id = format!("{}", unit.id);
            let record = fasta::Record::with_data(&id, &None, unit.seq.as_bytes());
            wtr.write_record(&record)?;
        }
        let mut reads = c_dir.clone();
        reads.push("reads.fa");
        let mut wtr = fasta::Writer::new(std::fs::File::create(&reads)?);
        for read in ds.raw_reads.iter() {
            let id = format!("{}", read.name);
            let record = fasta::Record::with_data(&id, &None, read.seq.as_bytes());
            wtr.write_record(&record)?;
        }
        let reference = reference.into_os_string().into_string().unwrap();
        let reads = reads.into_os_string().into_string().unwrap();
        (reference, reads)
    };
    let db_name = {
        let mut temp = c_dir.clone();
        temp.push("reference");
        temp.into_os_string().into_string().unwrap()
    };
    // Create database - train - align
    let lastdb = std::process::Command::new("lastdb")
        .args(&["-R", "00", "-Q", "0", &db_name, &reference])
        .output()?;
    if !lastdb.status.success() {
        panic!("lastdb-{}", String::from_utf8_lossy(&lastdb.stderr));
    }
    let p = format!("{}", p);
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
            "-f", "tab", "-P", &p, "-R", "00", "-Q", "0", "-p", &param, &db_name, &reads,
        ])
        .output()
        .unwrap();
    if !lastal.status.success() {
        panic!("lastal-{:?}", String::from_utf8_lossy(&lastal.stderr));
    }
    let alignments: Vec<_> = String::from_utf8_lossy(&lastal.stdout)
        .lines()
        .filter(|e| !e.starts_with("#"))
        .filter_map(|e| LastTAB::from_line(&e))
        .collect();
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    Ok(alignments)
}

fn distribute<'a>(alignments: &'a [LastTAB]) -> HashMap<String, Vec<&'a LastTAB>> {
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for alignment in alignments {
        let q_name = alignment.seq2_name().to_string();
        buckets.entry(q_name).or_default().push(alignment);
    }
    buckets
}

use definitions::{Edge, EncodedRead, Node, Op, RawRead, Unit};
fn encode(read: &RawRead, alignments: &[&LastTAB], units: &[Unit]) -> Option<EncodedRead> {
    let mut nodes: Vec<_> = alignments
        .iter()
        .filter(|aln| aln.seq1_matchlen() > aln.seq1_len() * 98 / 100)
        .filter_map(|aln| encode_alignment(aln, units, read))
        .collect();
    nodes.sort_by_key(|e| e.position_from_start);
    let edges: Vec<_> = nodes
        .windows(2)
        .map(|w| Edge::from_nodes(w, read.seq()))
        .collect();
    let leading_gap = nodes.first()?.position_from_start;
    let trailing_gap = {
        let last = nodes.last()?;
        read.seq().len() - last.position_from_start + consumed_reference_length(&last.cigar)
    };
    Some(EncodedRead {
        original_length: read.seq().len(),
        id: read.id,
        edges,
        nodes,
        leading_gap,
        trailing_gap,
    })
}

fn encode_alignment(aln: &LastTAB, _units: &[Unit], read: &RawRead) -> Option<Node> {
    let position = aln.seq2_start_from_forward();
    let unit: u64 = aln.seq1_name().parse().ok()?;
    let is_forward = aln.seq2_direction().is_forward();
    let seq = {
        let start = aln.seq2_start();
        let end = start + aln.seq2_matchlen();
        let seq = if is_forward {
            read.seq().to_vec()
        } else {
            bio_utils::revcmp(read.seq())
        };
        seq[start..end].to_vec()
    };
    let cigar = convert_aln_to_cigar(aln);
    Some(Node {
        position_from_start: position,
        unit,
        cluster: 0,
        seq: String::from_utf8_lossy(&seq).to_string(),
        is_forward,
        cigar,
    })
}

fn convert_aln_to_cigar(aln: &lasttab::LastTAB) -> Vec<Op> {
    let mut cigar = if aln.seq1_start_from_forward() != 0 {
        vec![Op::Del(aln.seq1_start_from_forward())]
    } else {
        vec![]
    };
    cigar.extend(aln.alignment().into_iter().map(|op| match op {
        lasttab::Op::Seq1In(l) => Op::Ins(l),
        lasttab::Op::Seq2In(l) => Op::Del(l),
        lasttab::Op::Match(l) => Op::Match(l),
    }));
    let reflen = consumed_reference_length(&cigar);
    assert!(reflen <= aln.seq1_len(), "{} > {}", reflen, aln.seq1_len());
    if aln.seq1_len() > reflen {
        cigar.push(Op::Del(aln.seq1_len() - reflen))
    }
    cigar
}

fn consumed_reference_length(cigar: &[Op]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            Op::Match(l) | Op::Del(l) => *l,
            Op::Ins(_) => 0,
        })
        .sum::<usize>()
}

use bio_utils::lasttab;
use bio_utils::lasttab::LastTAB;
use rayon::prelude::*;
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
            .par_iter()
            .filter_map(|read| {
                let alns = alignments_each_reads.get(&read.name)?;
                encode(read, alns, &self.selected_chunks)
            })
            .collect();
        debug!("Encoded {} reads.", encoded_reads.len());
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
            let id = read.name.to_string();
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
        .filter(|e| !e.starts_with('#'))
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
pub fn encode(read: &RawRead, alignments: &[&LastTAB], units: &[Unit]) -> Option<EncodedRead> {
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for &aln in alignments {
        let r_name = aln.seq1_name().to_string();
        let q_direction = aln.seq2_direction().is_forward();
        buckets.entry((r_name, q_direction)).or_default().push(aln);
    }
    let mut buckets: Vec<_> = buckets.values().collect();
    buckets.sort_by_key(|alns| alns.iter().map(|aln| aln.seq1_start()).min().unwrap_or(0));
    let mut nodes: Vec<_> = buckets
        .iter()
        .filter_map(|alns| encode_alignment(alns, units, read))
        .collect();
    nodes.sort_by_key(|e| e.position_from_start);
    let edges: Vec<_> = nodes
        .windows(2)
        .map(|w| Edge::from_nodes(w, read.seq()))
        .collect();
    let leading_gap = nodes.first()?.position_from_start;
    let trailing_gap = {
        let last = nodes.last()?;
        read.seq().len() - last.position_from_start - last.seq.len()
    };
    let len = nodes.iter().map(|n| n.seq.len()).sum::<usize>() as i64;
    let edge_len = edges.iter().map(|n| n.offset).sum::<i64>();
    let chunked_len = (len + edge_len) as usize;
    assert_eq!(read.seq().len(), chunked_len + leading_gap + trailing_gap);
    Some(EncodedRead {
        original_length: read.seq().len(),
        id: read.id,
        edges,
        nodes,
        leading_gap,
        trailing_gap,
    })
}

fn encode_alignment(alns: &[&LastTAB], units: &[Unit], read: &RawRead) -> Option<Node> {
    let aln = alns.get(0)?;
    let is_forward = aln.seq2_direction().is_forward();
    let seq = if is_forward {
        read.seq().to_vec()
    } else {
        bio_utils::revcmp(read.seq())
    };
    let unit_id: u64 = aln.seq1_name().parse().ok()?;
    let unit = units.iter().find(|u| u.id == unit_id)?;
    encode_alignment_by_chaining(alns, unit, &seq).map(|(position_from_start, seq, cigar)| {
        if is_forward {
            Node {
                position_from_start,
                unit: unit_id,
                cluster: 0,
                seq,
                is_forward,
                cigar,
            }
        } else {
            let position_from_start = read.seq().len() - position_from_start - seq.len();
            Node {
                position_from_start,
                unit: unit_id,
                cluster: 0,
                seq,
                is_forward,
                cigar,
            }
        }
    })
}

#[derive(Debug, Clone)]
enum DAGNode<'a> {
    Aln(&'a LastTAB),
    Start(usize),
    End(usize),
}

// Query start, Query sequence, Cigar.
type Alignment = (usize, String, Vec<Op>);
fn encode_alignment_by_chaining(alns: &[&LastTAB], unit: &Unit, read: &[u8]) -> Option<Alignment> {
    assert!(!alns.is_empty());
    // Compute anchor position.
    let (query_start, query_end) = get_read_range(alns, unit, read)?;
    if query_end - query_start < unit.seq().len() * 8 / 10 {
        return None;
    }
    let mut dag_nodes = vec![DAGNode::Start(query_start), DAGNode::End(query_end)];
    dag_nodes.extend(alns.iter().map(|&a| DAGNode::Aln(a)));
    let chain = chaining(&dag_nodes, unit.seq().len());
    // Check if this chain covers sufficient fraction of the unit.
    let cover_length = chain
        .iter()
        .map(|chain| match chain {
            DAGNode::Aln(x) => x.seq1_matchlen(),
            _ => 0,
        })
        .sum::<usize>();
    if cover_length < unit.seq().len() * 9 / 10 {
        return None;
    }
    let (mut q_pos, mut r_pos) = (query_start, 0);
    let mut cigar = vec![];
    match *chain[0] {
        DAGNode::Start(x) => assert_eq!(x, q_pos),
        _ => panic!(),
    }
    for node in chain.iter().skip(1) {
        let (q_target, r_target) = match node {
            &&DAGNode::End(x) => (x, unit.seq().len()),
            DAGNode::Aln(x) => (x.seq2_start(), x.seq1_start()),
            _ => panic!(),
        };
        // Make them as parameters
        let (query, refr) = (&read[q_pos..q_target], &unit.seq()[r_pos..r_target]);
        cigar.extend(alignment(query, refr, 1, -1, -1));
        if let DAGNode::Aln(x) = node {
            cigar.extend(convert_aln_to_cigar(x));
            r_pos = x.seq1_start() + x.seq1_matchlen();
            q_pos = x.seq2_start() + x.seq2_matchlen();
        }
    }
    cigar.reverse();
    let mut query_start = query_start;
    while let Some(&Op::Ins(l)) = cigar.last() {
        query_start += l;
        cigar.pop();
    }
    cigar.reverse();
    let mut query_end = query_end;
    while let Some(&Op::Ins(l)) = cigar.last() {
        query_end -= l;
        cigar.pop();
    }
    // if query_end - query_start > unit.seq().len() * 11 / 10 {
    //     debug!("LEN:{}", query_end - query_start);
    //     for aln in alns.iter() {
    //         debug!("DUMP\t{}", aln);
    //     }
    //     let query = &read[query_start..query_end];
    //     let (q, op, r) = recover(query, unit.seq(), &cigar);
    //     for ((q, op), r) in q.chunks(100).zip(op.chunks(100)).zip(r.chunks(100)) {
    //         eprintln!("{}", String::from_utf8_lossy(q));
    //         eprintln!("{}", String::from_utf8_lossy(op));
    //         eprintln!("{}", String::from_utf8_lossy(r));
    //     }
    // }
    assert_eq!(consumed_reference_length(&cigar), unit.seq().len());
    let query = String::from_utf8_lossy(&read[query_start..query_end]).to_string();
    Some((query_start, query, cigar))
}

pub fn recover(query: &[u8], refr: &[u8], ops: &[Op]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let (mut q, mut al, mut r) = (vec![], vec![], vec![]);
    let (mut q_pos, mut r_pos) = (0, 0);
    for op in ops {
        match *op {
            Op::Match(l) => {
                al.extend(
                    query[q_pos..q_pos + l]
                        .iter()
                        .zip(&refr[r_pos..r_pos + l])
                        .map(|(x, y)| if x == y { b'|' } else { b'X' }),
                );
                q.extend(query[q_pos..q_pos + l].iter().copied());
                r.extend(refr[r_pos..r_pos + l].iter().copied());
                q_pos += l;
                r_pos += l;
            }
            Op::Del(l) => {
                al.extend(vec![b' '; l]);
                q.extend(vec![b' '; l]);
                r.extend(refr[r_pos..r_pos + l].iter().copied());
                r_pos += l;
            }
            Op::Ins(l) => {
                al.extend(vec![b' '; l]);
                q.extend(query[q_pos..q_pos + l].iter().copied());
                r.extend(vec![b' '; l]);
                q_pos += l;
            }
        }
    }
    (q, al, r)
}

fn get_read_range(alns: &[&LastTAB], unit: &Unit, read: &[u8]) -> Option<(usize, usize)> {
    let query_start = {
        let first = alns.iter().min_by_key(|aln| aln.seq1_start())?;
        let remaining = first.seq1_start();
        let query_position = first.seq2_start();
        query_position.max(2 * remaining) - 2 * remaining
    };
    let query_end = {
        let last = alns
            .iter()
            .max_by_key(|aln| aln.seq1_start() + aln.seq1_matchlen())?;
        let remaining = unit.seq().len() - (last.seq1_start() + last.seq1_matchlen());
        let query_position = last.seq2_start() + last.seq2_matchlen();
        (query_position + 2 * remaining).min(read.len())
    };
    Some((query_start, query_end))
}

fn chaining<'a, 'b>(nodes: &'b [DAGNode<'a>], unit_len: usize) -> Vec<&'b DAGNode<'a>> {
    let edges: Vec<Vec<(usize, i64)>> = nodes
        .iter()
        .map(|u| {
            nodes
                .iter()
                .enumerate()
                .filter_map(|(idx, v)| compute_edge(u, v, unit_len).map(|w| (idx, w)))
                .collect()
        })
        .collect();
    let order = topological_sort(&edges);
    let (mut dp, mut parent) = (vec![0; nodes.len()], vec![0; nodes.len()]);
    for i in order {
        for &(j, score) in edges[i].iter() {
            let aln_j = match nodes[j] {
                DAGNode::Aln(aln) => aln.score() as i64,
                _ => 0,
            };
            let from_i_to_j = dp[i] + score + aln_j;
            if dp[j] <= from_i_to_j {
                dp[j] = from_i_to_j;
                parent[j] = i;
            }
        }
    }
    // i <= 1, tracing back from the end node.
    assert!(match nodes[1] {
        DAGNode::End(_) => true,
        _ => false,
    });
    let mut path = vec![];
    let mut current_node = 1;
    while current_node != 0 {
        path.push(&nodes[current_node]);
        current_node = parent[current_node];
    }
    path.push(&nodes[current_node]);
    path.reverse();
    path
}

// Compute edge score between from u to v. If no edge possible, return None.
fn compute_edge<'a>(u: &DAGNode<'a>, v: &DAGNode<'a>, unit_len: usize) -> Option<i64> {
    match u {
        DAGNode::End(_) => None,
        DAGNode::Start(x) => match v {
            DAGNode::Aln(aln) => Some(-((aln.seq2_start() - x + aln.seq1_start()) as i64)),
            _ => None,
        },
        DAGNode::Aln(aln1) => match v {
            DAGNode::End(x) => {
                let refr_remain = unit_len - aln1.seq1_start() - aln1.seq1_matchlen();
                let query_remain = x - aln1.seq2_start() - aln1.seq2_matchlen();
                Some(-((refr_remain + query_remain) as i64))
            }
            DAGNode::Start(_) => None,
            DAGNode::Aln(aln2) => {
                let u_refr_end = aln1.seq1_start() + aln1.seq1_matchlen();
                let u_query_end = aln1.seq2_start() + aln1.seq2_matchlen();
                let v_refr_start = aln2.seq1_start();
                let v_query_start = aln2.seq2_start();
                if u_refr_end <= v_refr_start && u_query_end <= v_query_start {
                    Some(-((v_refr_start - u_refr_end + v_query_start - u_query_end) as i64))
                } else {
                    None
                }
            }
        },
    }
}

// Return nodes in topological order. Note that the graph IS connected.
fn topological_sort(edges: &[Vec<(usize, i64)>]) -> Vec<usize> {
    let len = edges.len();
    let mut arrived = vec![false; len];
    let mut stack = vec![0];
    let mut order = vec![];
    'dfs: while !stack.is_empty() {
        let node = *stack.last().unwrap();
        if !arrived[node] {
            arrived[node] = true;
        }
        for &(to, _) in &edges[node] {
            if !arrived[to] {
                stack.push(to);
                continue 'dfs;
            }
        }
        let last = stack.pop().unwrap();
        order.push(last);
    }
    order.reverse();
    order
}

fn alignment(query: &[u8], refr: &[u8], mat: i64, mism: i64, gap: i64) -> Vec<Op> {
    let mut dp = vec![vec![0; query.len() + 1]; refr.len() + 1];
    for j in 0..=query.len() {
        dp[0][j] = j as i64 * gap;
    }
    for (i, row) in dp.iter_mut().enumerate() {
        row[0] = i as i64 * gap;
        // for i in 0..=refr.len() {
        // dp[i][0] = i as i64 * gap;
    }
    for (i, r) in refr.iter().enumerate() {
        for (j, q) in query.iter().enumerate() {
            let match_score = if r == q { mat } else { mism };
            let max = (dp[i][j] + match_score)
                .max(dp[i][j + 1] + gap)
                .max(dp[i + 1][j] + gap);
            dp[i + 1][j + 1] = max;
        }
    }
    // Traceback.
    let (mut q_pos, mut r_pos) = (query.len(), refr.len());
    let mut ops = vec![];
    while 0 < q_pos && 0 < r_pos {
        let current = dp[r_pos][q_pos];
        if current == dp[r_pos - 1][q_pos] + gap {
            ops.push(2);
            r_pos -= 1;
        } else if current == dp[r_pos][q_pos - 1] + gap {
            ops.push(1);
            q_pos -= 1;
        } else {
            let match_score = if query[q_pos - 1] == refr[r_pos - 1] {
                mat
            } else {
                mism
            };
            assert_eq!(current, dp[r_pos - 1][q_pos - 1] + match_score);
            ops.push(0);
            q_pos -= 1;
            r_pos -= 1;
        }
    }
    for _ in 0..q_pos {
        ops.push(1);
    }
    for _ in 0..r_pos {
        ops.push(2);
    }
    compress(ops)
}

fn compress(mut ops: Vec<u8>) -> Vec<Op> {
    let mut cigar = vec![];
    while !ops.is_empty() {
        let last = ops.pop().unwrap();
        let mut count = 1;
        while let Some(&res) = ops.last() {
            if res == last {
                count += 1;
                ops.pop();
            } else {
                break;
            }
        }
        match last {
            0 => cigar.push(Op::Match(count)),
            1 => cigar.push(Op::Ins(count)),
            2 => cigar.push(Op::Del(count)),
            _ => panic!(),
        }
    }
    cigar
}

fn convert_aln_to_cigar(aln: &lasttab::LastTAB) -> Vec<Op> {
    aln.alignment()
        .into_iter()
        .map(|op| match op {
            lasttab::Op::Seq1In(l) => Op::Ins(l),
            lasttab::Op::Seq2In(l) => Op::Del(l),
            lasttab::Op::Match(l) => Op::Match(l),
        })
        .collect()
    // let mut cigar = if aln.seq1_start_from_forward() != 0 {
    //     vec![Op::Del(aln.seq1_start_from_forward())]
    // } else {
    //     vec![]
    // };
    // cigar.extend(aln.alignment().into_iter().map(|op| match op {
    //     lasttab::Op::Seq1In(l) => Op::Ins(l),
    //     lasttab::Op::Seq2In(l) => Op::Del(l),
    //     lasttab::Op::Match(l) => Op::Match(l),
    // }));
    // let reflen = consumed_reference_length(&cigar);
    // assert!(reflen <= aln.seq1_len(), "{} > {}", reflen, aln.seq1_len());
    // if aln.seq1_len() > reflen {
    //     cigar.push(Op::Del(aln.seq1_len() - reflen))
    // }
    // cigar
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn works() {}
    #[test]
    fn alignment_check() {
        let query = b"AAAAA";
        let reference = b"AAAAA";
        let res = alignment(query, reference, 1, -1, -1);
        assert_eq!(res, vec![Op::Match(query.len())]);
        let query = b"AAAAA";
        let reference = b"AACAA";
        let res = alignment(query, reference, 1, -1, -1);
        assert_eq!(res, vec![Op::Match(query.len())]);
        let query = b"AACCAAAA";
        let refer = b"AAGCAA";
        let res = alignment(query, refer, 1, -1, -1);
        assert_eq!(res, vec![Op::Match(refer.len()), Op::Ins(2)]);
        let query = b"ACGCGCGCAA";
        let refer = b"GCGCGC";
        let res = alignment(query, refer, 1, -1, -1);
        assert_eq!(res, vec![Op::Ins(2), Op::Match(6), Op::Ins(2)]);
        let query = b"GCGCGC";
        let refer = b"ACGCGCGCAA";
        let res = alignment(query, refer, 1, -1, -1);
        assert_eq!(res, vec![Op::Del(2), Op::Match(6), Op::Del(2)]);
        let query = b"CGCTGCGCAAAAA";
        let refer = b"AAAAAGCGCGCT";
        let res = alignment(query, refer, 1, -1, -1);
        let ans = vec![
            Op::Match(1),
            Op::Del(4),
            Op::Match(2),
            Op::Ins(1),
            Op::Match(5),
            Op::Ins(4),
        ];
        assert_eq!(res, ans)
    }
}

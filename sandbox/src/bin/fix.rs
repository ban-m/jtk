use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for read in ds.encoded_reads.iter_mut() {
        let mut nodes = Vec::with_capacity(read.nodes.len());
        nodes.append(&mut read.nodes);
        use haplotyper::encode::{
            nodes_to_encoded_read, remove_overlapping_encoding, remove_slippy_alignment,
        };
        let mut seq = read.recover_raw_read();
        seq.iter_mut().for_each(u8::make_ascii_uppercase);
        nodes.sort_by_key(|n| (n.chunk, n.position_from_start));
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        nodes = remove_overlapping_encoding(nodes);
        nodes = remove_contained_encoding(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, &seq).unwrap();
    }
    println!("{}", serde_json::ser::to_string_pretty(&ds).unwrap());
    Ok(())
}

fn remove_contained_encoding(mut nodes: Vec<Node>) -> Vec<Node> {
    'outer: loop {
        let len = nodes.len();
        for i in 0..len.max(1) - 1 {
            let p_start = nodes[i].position_from_start;
            let p_end = nodes[i].query_length() + p_start;
            let a_start = nodes[i + 1].position_from_start;
            let a_end = nodes[i + 1].query_length() + a_start;
            let p_contains_a = p_start <= a_start && a_end < p_end;
            let a_contains_p = a_start <= p_start && p_end < a_end;
            if p_contains_a || a_contains_p {
                let p_diff: usize = nodes[i]
                    .cigar
                    .iter()
                    .map(|op| match *op {
                        Op::Match(_) => 0,
                        Op::Del(l) => l,
                        Op::Ins(l) => l,
                    })
                    .sum();
                let a_diff: usize = nodes[i + 1]
                    .cigar
                    .iter()
                    .map(|op| match *op {
                        Op::Match(_) => 0,
                        Op::Del(l) => l,
                        Op::Ins(l) => l,
                    })
                    .sum();
                let p_len = nodes[i].query_length();
                let a_len = nodes[i + 1].query_length();
                let p_err = (p_diff as f64) / (p_len as f64);
                let a_err = (a_diff as f64) / (a_len as f64);
                let p_chunk = nodes[i].chunk;
                let a_chunk = nodes[i + 1].chunk;
                eprintln!("{p_start}-{p_end}({p_err:.3}){p_chunk}   {a_start}-{a_end}({a_err:.3}{a_chunk})");
                if p_err < a_err {
                    nodes.remove(i + 1);
                } else {
                    nodes.remove(i);
                }
                continue 'outer;
            }
        }
        break;
    }
    nodes
}

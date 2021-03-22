use definitions::*;
use serde_json;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    eprintln!("{:?}", std::time::Instant::now() - start);
    let mut degree: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|u| (u.id, HashSet::new()))
        .collect();
    for read in ds.encoded_reads.iter() {
        for edge in read.edges.iter() {
            degree.get_mut(&edge.from).unwrap().insert(edge.to);
            degree.get_mut(&edge.to).unwrap().insert(edge.from);
        }
    }
    let removed_unit: HashSet<_> = degree
        .iter()
        .filter_map(|(unit, degree)| (degree.len() > 10).then(|| *unit))
        .collect();
    eprintln!("Remove:{:?}", removed_unit);
    for read in ds.encoded_reads.iter_mut() {
        while let Some(idx) = read
            .nodes
            .iter()
            .position(|n| removed_unit.contains(&n.unit))
        {
            remove_ith_node(read, idx);
        }
    }
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| definitions::Assignment::new(r.id, 0))
        .collect();
    use haplotyper::Assemble;
    let assemble_config = haplotyper::AssembleConfig::new(1, 100, false);
    eprintln!("Start assembling {} reads", ds.encoded_reads.len());
    eprintln!("Assembled reads.");
    let gfa = ds.assemble_as_gfa(&assemble_config);
    println!("{}", gfa);
    Ok(())
}

pub fn remove_ith_node(read: &mut EncodedRead, i: usize) {
    assert!(i < read.nodes.len());
    let node = read.nodes.remove(i);
    let seq = if node.is_forward {
        node.seq().to_vec()
    } else {
        bio_utils::revcmp(node.seq())
    };
    if read.nodes.is_empty() {
        assert!(read.edges.is_empty());
        read.leading_gap.extend(seq);
        read.leading_gap.extend_from_slice(&read.trailing_gap);
        read.trailing_gap.clear();
        return;
    }
    if i == 0 {
        let edge = read.edges.remove(i);
        read.leading_gap.extend(seq);
        for _ in 0..(-edge.offset).max(0) {
            read.leading_gap.pop();
        }
        read.leading_gap.extend_from_slice(edge.label.as_bytes());
    } else if i == read.nodes.len() {
        let edge = read.edges.remove(i - 1);
        let mut trailing = edge.label.as_bytes().to_vec();
        trailing.extend_from_slice(&seq[(-edge.offset).max(0) as usize..]);
        trailing.extend_from_slice(&read.trailing_gap);
        read.trailing_gap = trailing;
    } else {
        let edge = read.edges.remove(i);
        read.edges[i - 1].to = edge.to;
        let mut new_label = read.edges[i - 1].label.as_bytes().to_vec();
        new_label.extend_from_slice(&seq[(-read.edges[i - 1].offset).max(0) as usize..]);
        for _ in 0..(-edge.offset).max(0) {
            new_label.pop();
        }
        new_label.extend_from_slice(edge.label.as_bytes());
        read.edges[i - 1].label = String::from_utf8(new_label).unwrap();
        read.edges[i - 1].offset += seq.len() as i64 + edge.offset;
        assert_eq!(
            read.edges[i - 1].offset.max(0),
            read.edges[i - 1].label.as_bytes().len() as i64
        );
    }
}

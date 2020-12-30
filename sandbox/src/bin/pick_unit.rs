// use bio_utils::fasta::Record;
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let mut unit_count: HashMap<_, u32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *unit_count.entry(node.unit).or_default() += 1;
        }
    }
    let mut unit_count: Vec<_> = unit_count.iter().collect();
    unit_count.sort_by_key(|x| x.0);
    println!("ID\tCount");
    for (id, count) in unit_count.iter() {
        println!("{}\t{}", id, count);
    }
    // let unit: u64 = args[2].parse().unwrap();
    // let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // let reads: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .filter_map(|read| {
    //         let name = id2name[&read.id];
    //         read.nodes
    //             .iter()
    //             .filter(|c| c.unit == unit)
    //             .next()
    //             .map(|n| Record::with_data(&name, &None, n.seq.as_bytes()))
    //     })
    //     .collect();
    // for read in reads {
    //     println!("{}", read);
    // }
}

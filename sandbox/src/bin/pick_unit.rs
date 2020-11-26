use bio_utils::fasta::Record;
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let unit: u64 = args[2].parse().unwrap();
    let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter_map(|read| {
            let name = id2name[&read.id];
            read.nodes
                .iter()
                .filter(|c| c.unit == unit)
                .next()
                .map(|n| Record::with_data(&name, &None, n.seq.as_bytes()))
        })
        .collect();
    for read in reads {
        println!("{}", read);
    }
}

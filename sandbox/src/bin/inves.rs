#![allow(unused_imports)]
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let (read_error_rate, unit_error_rate, sd) =
        haplotyper::encode::deletion_fill::estimate_error_rate(&ds, 0.15);
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    eprintln!("{}", sd);
    // let mut ok = 0;
    // println!("read\tunit\terror\texpected\tthr");
    // for read in ds.encoded_reads.iter() {
    //     for node in read.nodes.iter() {
    //         let (_, aln, _) = node.recover(units[&node.unit]);
    //         let diff = aln.iter().filter(|&&x| x != b'|').count();
    //         let error_rate = diff as f64 / aln.len() as f64;
    //         let expected = read_error_rate[read.id as usize] + unit_error_rate[node.unit as usize];
    //         let threshold = expected + 10f64 * sd;
    //         println!(
    //             "{}\t{}\t{}\t{}\t{}",
    //             read.id, node.unit, error_rate, expected, threshold
    //         );
    //         ok += (error_rate < threshold) as usize;
    //     }
    // }
    // let total: usize = ds.encoded_reads.iter().map(|r| r.nodes.len()).sum();
    // eprintln!("{}\t{}", ok, total);
    for read in ds.encoded_reads.iter() {
        let error = read_error_rate[read.id as usize];
        let errors: f64 = read
            .nodes
            .iter()
            .map(|node| {
                let (_, aln, _) = node.recover(units[&node.unit]);
                let diff = aln.iter().filter(|&&x| x != b'|').count();
                diff as f64 / aln.len() as f64
            })
            .sum();
        let ave = errors / read.nodes.len() as f64;
        let len = read.nodes.len();
        let orig = read.original_length;
        println!("{}\t{}\t{}\t{}\t{}\tRead", read.id, error, ave, len, orig);
    }
    for chunk in ds.selected_chunks.iter() {
        let error = unit_error_rate[chunk.id as usize];
        println!("{}\t{}\tUnit", chunk.id, error);
    }
    Ok(())
}

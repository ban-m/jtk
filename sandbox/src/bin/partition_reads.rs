use bio_utils::fastq;
use std::collections::HashMap;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::BufWriter;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let partitions: HashMap<_, _> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .map_while(Result::ok)
        .skip(1)
        .filter_map(|line| {
            let mut line = line.split_whitespace();
            let readid = line.next()?.to_string();
            let hap = match line.next()? {
                "none" => 0,
                "H1" => 1,
                "H2" => 2,
                _ => panic!(),
            };
            Some((readid, hap))
        })
        .collect();
    let reads = fastq::parse_into_vec(&args[2])?;
    let mut hap1_wtr = std::fs::File::create(&args[3]).map(BufWriter::new)?;
    let mut hap2_wtr = std::fs::File::create(&args[4]).map(BufWriter::new)?;
    for read in reads.iter() {
        match partitions.get(read.id()) {
            Some(0) | None => {
                writeln!(hap1_wtr, "{read}")?;
                writeln!(hap2_wtr, "{read}")?
            }
            Some(1) => writeln!(hap1_wtr, "{read}")?,
            Some(2) => writeln!(hap2_wtr, "{read}")?,
            _ => panic!(),
        }
    }
    Ok(())
}

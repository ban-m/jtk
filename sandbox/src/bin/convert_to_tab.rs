use bio_utils::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::*;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut header = HashMap::new();
    let records: Vec<_> = File::open(&args[1])
        .map(BufReader::new)
        .unwrap()
        .lines()
        .filter_map(|l| l.ok())
        .filter_map(|line| {
            if line.starts_with("@SQ") {
                let mut line = line.split('\t');
                line.next();
                let name = line.next()?.split(':').nth(1)?.to_string();
                let length: usize = line.next()?.split(':').nth(1)?.parse().ok()?;
                header.insert(name, length);
                None
            } else {
                sam::Sam::new(&line)
            }
        })
        .collect();
    let reference = bio_utils::fasta::parse_into_vec(&args[2]).unwrap();
    let reads = bio_utils::fasta::parse_into_vec(&args[3]).unwrap();
    for (k, v) in header.iter().take(10) {
        eprintln!("{}\t{}", k, v);
    }
    let mut result: HashMap<_, Vec<_>> = HashMap::new();
    for record in records {
        match lasttab::try_from(&record, &header) {
            Ok(tab) => result
                .entry(tab.seq1_name().to_string())
                .or_default()
                .push(tab),
            Err(why) => eprintln!("{}", why),
        }
    }
    for (refname, alns) in result {
        if let Some(template) = reference.iter().find(|r| r.id() == refname) {
            let pileup = haplotyper::assemble::pileup::Pileups::convert_into_pileup_raw(
                &alns, template, &reads,
            );
            pileup.dump();
        }
    }
}

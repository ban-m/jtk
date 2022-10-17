const CHUNK_SIZE: usize = 1_000;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::collections::HashMap;
    let target_regions: HashMap<_, _> = args[1]
        .split_whitespace()
        .map(|region| {
            if region.contains(':') {
                let (region, range) = region.rsplit_once(':').unwrap();
                let (start, end) = range.rsplit_once('-').unwrap();
                let start: usize = start.parse().unwrap();
                let end: usize = end.parse().unwrap();
                (region.to_string(), Some((start, end)))
            } else {
                (region.to_string(), None)
            }
        })
        .collect();
    for (name, range) in target_regions.iter() {
        eprintln!("{name},{range:?}");
    }
    let stdin = std::io::stdin();
    use std::io::*;
    let mut coverages: HashMap<_, Vec<_>> = HashMap::new();
    for line in std::io::BufReader::new(stdin.lock())
        .lines()
        .filter_map(|l| l.ok())
    {
        let mut line = line.split_whitespace();
        let ctgname = line.next().unwrap();
        let position: usize = line.next().unwrap().parse().unwrap();
        let cov: usize = line.next().unwrap().parse().unwrap();
        if target_regions.contains_key(ctgname) {
            let chunk_id = position / CHUNK_SIZE;
            let slot = coverages.entry(ctgname.to_string()).or_default();
            if slot.len() <= chunk_id {
                let len = (chunk_id + 1) - slot.len();
                slot.extend(std::iter::repeat((0, 0)).take(len));
            }
            slot[chunk_id].0 += cov;
            slot[chunk_id].1 += 1;
        }
    }
    println!("Contig\tPosition\tCoverage");
    for (ctgname, coverages) in coverages.iter() {
        let region = target_regions[ctgname];
        let coverages = coverages
            .iter()
            .enumerate()
            .map(|(i, x)| (i * CHUNK_SIZE, x));
        if let Some((start, end)) = region {
            for (pos, (total, num)) in coverages.filter(|&(pos, _)| start <= pos && pos <= end) {
                println!("{ctgname}\t{}\t{}", pos - start, total / num);
            }
        } else {
            for (pos, (total, num)) in coverages {
                println!("{ctgname}\t{pos}\t{}", total / num);
            }
        }
    }
    Ok(())
}

use definitions::*;
use std::collections::HashSet;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let tigs = std::fs::File::open(&args[2]).map(BufReader::new)?;
    use std::collections::HashMap;
    use std::io::BufRead;
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let tigs: Vec<_> = tigs
        .lines()
        .filter_map(|x| x.ok())
        .map(|line| {
            let mut line = line.split('\t');
            let id = line.next().unwrap().to_string();
            let copy_num: usize = line.next().unwrap().parse().unwrap();
            let clusters: HashSet<_> = line
                .filter_map(|x| {
                    let mut elm = x.split('-');
                    let unit: u64 = elm.next()?.parse().unwrap();
                    let cluster: u64 = elm.next()?.parse().unwrap();
                    Some((unit, cluster))
                })
                .collect();
            (id, copy_num, clusters)
        })
        .collect();
    use sandbox::IS_MOCK;
    let summaries: Vec<_> = tigs.iter().map(|x| x.2.clone()).collect();
    let mut assign_contigs = vec![[0; 2]; summaries.len()];
    use rayon::prelude::*;
    let dists: Vec<_> = ds
        .encoded_reads
        .par_iter()
        .map(|read| {
            let nodes: Vec<_> = read.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
            let dist = haplotyper::assemble::distribute(&nodes, &summaries);
            (read, dist)
        })
        .collect();
    for (read, dist) in dists {
        let answer = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        let mut counts = vec![0; summaries.len()];
        for d in dist {
            counts[d] += 1;
        }
        let (asn, _) = counts.iter().enumerate().max_by_key(|x| x.1).unwrap();
        assign_contigs[asn][answer] += 1;
    }
    for ((id, copy_num, _), counts) in tigs.iter().zip(assign_contigs) {
        println!("{id}\t{copy_num}\t{}\t{}", counts[0], counts[1]);
    }
    Ok(())
}

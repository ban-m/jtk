use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() {
    let start = std::time::Instant::now();
    let ds = get_input_file().unwrap();
    let end = std::time::Instant::now();
    eprintln!("{:?}", end - start);
    let k = 4;
    let mut count: HashMap<_, HashMap<_, usize>> = HashMap::new();
    let is_hap_a: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| {
            let desc = if read.desc.contains("252v2") { 1 } else { 0 };
            (read.id, desc)
        })
        .collect();
    for read in ds.encoded_reads.iter().filter(|r| r.nodes.len() >= 4) {
        let cluster = is_hap_a[&read.id];
        for node in read.nodes.iter() {
            *count
                .entry(node.unit)
                .or_default()
                .entry((cluster, node.cluster))
                .or_default() += 1;
        }
    }
    let cluster_num: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|unit| (unit.id, unit.cluster_num))
        .collect();
    let mut count: Vec<_> = count.into_iter().collect();
    count.sort_by_key(|x| x.0);
    for (slot, result) in count {
        // let mut line = vec![];
        let cluster_num = cluster_num[&slot] as u64;
        println!("UNIT\t{}\t{}", slot, cluster_num);
        for ans in 0..cluster_num {
            for pred in 0..cluster_num {
                let count = match result.get(&(ans, pred)) {
                    Some(count) => count,
                    None => &0,
                };
                print!("\t{}", count);
            }
            println!();
        }
        println!();
        // let total = result.values().sum::<usize>();
        // let line = line.join("\t");
        // println!("UNIT\t{}\t{}\t{}\t{}", slot, line, total, cluster_num);
    }
    let mut count = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for kmer in read.nodes.windows(k) {
            let mut kmer: Vec<_> = kmer.iter().map(|u| (u.unit, u.cluster)).collect();
            if kmer.last().unwrap() <= kmer.first().unwrap() {
                kmer.reverse();
            }
            *count.entry(kmer).or_default() += 1;
        }
    }
    let count: Vec<_> = count.values().copied().collect();
    let hist = histgram_viz::Histgram::new(&count);
    println!("{}", hist.format(20, 20));
    let mut kmer_in_hap_d: HashMap<Vec<(u64, u64)>, usize> = HashMap::new();
    let mut kmer_in_hap_a: HashMap<Vec<(u64, u64)>, usize> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let is_in_a = is_hap_a[&read.id] == 1;
        for kmer in read.nodes.windows(k) {
            let mut kmer: Vec<_> = kmer.iter().map(|n| (n.unit, n.cluster as u64)).collect();
            if kmer.last().unwrap() < kmer.first().unwrap() {
                kmer.reverse();
            }
            if is_in_a {
                *kmer_in_hap_a.entry(kmer).or_default() += 1;
            } else {
                *kmer_in_hap_d.entry(kmer).or_default() += 1;
            }
        }
    }
    println!("{}\t{}", kmer_in_hap_d.len(), kmer_in_hap_a.len());
    for (kmer, count_d) in kmer_in_hap_d.iter() {
        if let Some(count_a) = kmer_in_hap_a.get(kmer) {
            println!("{:?}\t{}\t{}", kmer, count_d, count_a);
        }
    }
}

fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

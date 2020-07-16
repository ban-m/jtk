use definitions::*;
use rayon::prelude::*;
#[macro_use]
extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("trace")).init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let units = dataset.selected_chunks.clone();
    use std::time::Instant;
    let start = Instant::now();
    let k: u32 = 7;
    let kmer_vec: Vec<_> = units
        .iter()
        .map(|unit| {
            let mut fingerprint = vec![false; 4usize.pow(k)];
            for kmer in unit.seq.as_bytes().windows(k as usize) {
                fingerprint[to_index(kmer)] = true;
            }
            fingerprint
        })
        .collect();
    println!("Finger:{:?}", Instant::now() - start);
    let result: Vec<_> = units
        .par_iter()
        .enumerate()
        .flat_map(|(idx, unit1)| {
            units
                .iter()
                .enumerate()
                .skip(idx + 1)
                .filter_map(|(jdx, unit2)| {
                    let share_k_mer = kmer_vec[idx]
                        .iter()
                        .zip(kmer_vec[jdx].iter())
                        .filter(|(&x, &y)| x && y)
                        .count();
                    let tot_k_mer = kmer_vec[idx]
                        .iter()
                        .zip(kmer_vec[jdx].iter())
                        .filter(|(&x, &y)| x || y)
                        .count();
                    let jaccard = share_k_mer as f64 / tot_k_mer as f64;
                    if jaccard > 0.7 {
                        let aln = alignment(unit1.seq.as_bytes(), unit2.seq.as_bytes());
                        Some((unit1.id, unit2.id, aln))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let end = Instant::now();
    println!("{}/{}", result.len(), units.len() * units.len());
    println!("{:?}", end - start);
    for (id1, id2, score) in result {
        println!("{}\t{}\t{}", id1, id2, score);
    }
    // use std::collections::HashMap;
    // let id_to_asn: HashMap<_, _> = dataset
    //     .assignments
    //     .iter()
    //     .map(|asn| (asn.id, asn.cluster))
    //     .collect();
    // let (unit, cluster, group) = (1517, 1, 1);
    // for read in dataset
    //     .encoded_reads
    //     .iter()
    //     .filter(|read| match id_to_asn.get(&read.id) {
    //         Some(&asn) if asn == group => read
    //             .nodes
    //             .iter()
    //             .any(|n| n.unit == unit && n.cluster == cluster),
    //         _ => false,
    //     })
    // {
    //     let units: Vec<_> = read
    //         .nodes
    //         .windows(2)
    //         .zip(read.edges.iter())
    //         .map(|(n, e)| format!("{}:{}-({})-", n[0].unit, n[0].cluster, e.offset))
    //         .collect();
    //     println!("{}", units.join(""));
    // }
    Ok(())
}

fn to_index(kmer: &[u8]) -> usize {
    kmer.iter().fold(0, |x, y| match y {
        b'A' => (x << 2),
        b'C' => (x << 2) + 1,
        b'G' => (x << 2) + 2,
        b'T' => (x << 2) + 3,
        _ => panic!(),
    })
}

fn alignment(seq1: &[u8], seq2: &[u8]) -> i32 {
    let mut dp = vec![vec![0; seq2.len() + 1]; seq1.len() + 1];
    for i in 1..seq1.len() + 1 {
        for j in 1..seq2.len() + 1 {
            let m = if seq1[i - 1] == seq2[j - 1] { 1 } else { -1 };
            dp[i][j] = (dp[i - 1][j] - 1)
                .max(dp[i][j - 1] - 1)
                .max(dp[i - 1][j - 1] + m);
        }
    }
    dp[seq1.len()][seq2.len()].max(0)
}

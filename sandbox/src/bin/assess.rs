use definitions::*;
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Write};
#[derive(Clone, Serialize, Deserialize)]
pub struct ERead {
    pub id: u64,
    pub path: Vec<Elm>,
    pub cluster: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Elm {
    pub unit: u64,
    pub cluster: usize,
}
use std::collections::HashMap;
/// Assess clustering goodness.
/// # Synopsis
/// ```bash
/// cargo run --release --bin assess -- ${JSON encoded file}
/// ```
/// would create local.tsv, global.tsv
fn main() -> std::io::Result<()> {
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let answer: HashMap<_, usize> = dataset
        .raw_reads
        .iter()
        .map(|r| {
            let cluster = if r.name.contains("hapA") {
                0
            } else if r.name.contains("hapB") {
                1
            } else if r.name.contains("hapC") {
                2
            } else {
                eprintln!("Unknown name:{}", r.name);
                3
            };
            (r.id, cluster)
        })
        .collect();
    local_assess(&answer, &dataset)?;
    global_assess(&answer, &dataset)?;
    Ok(())
}

fn local_assess(answer: &HashMap<u64, usize>, ds: &definitions::DataSet) -> std::io::Result<()> {
    let pileup = {
        let mut result: HashMap<_, Vec<_>> = HashMap::new();
        for eread in ds.encoded_reads.iter() {
            let ans = match answer.get(&eread.id) {
                Some(res) => res,
                None => continue,
            };
            for node in eread.nodes.iter() {
                result
                    .entry(node.unit)
                    .or_default()
                    .push((*ans, node.cluster));
            }
        }
        result
    };
    let mut out = BufWriter::new(std::fs::File::create("local.tsv")?);
    writeln!(&mut out, "Unit\tCluster\tPurity\tConsent")?;
    for (unit, predictions) in pileup {
        for (cluster, purity, consent) in unit_assess(&predictions) {
            let purity = if purity.is_nan() { 0. } else { purity };
            let consent = if consent.is_nan() { 0. } else { consent };
            writeln!(&mut out, "{}\t{}\t{}\t{}", unit, cluster, purity, consent)?;
        }
    }
    Ok(())
}

fn unit_assess(preds: &[(usize, u64)]) -> Vec<(usize, f64, f64)> {
    // There are three clusters, we know.
    let mut counts = [[0; 3]; 3];
    for &(ans, pred) in preds {
        counts[ans][pred as usize] += 1;
    }
    eprintln!();
    for cs in counts.iter() {
        let cs: Vec<_> = cs.iter().map(|e| format!("{}", e)).collect();
        eprintln!("{}", cs.join("\t"));
    }
    (0..3)
        .map(|cl| {
            let purity = {
                let tot = counts.iter().map(|x| x[cl]).sum::<usize>();
                let kth = counts.iter().map(|x| x[cl]).max().unwrap();
                kth as f64 / tot as f64
            };
            let consent = {
                let tot = counts[cl].iter().sum::<usize>();
                let kth = *counts[cl].iter().max().unwrap();
                kth as f64 / tot as f64
            };
            (cl, purity, consent)
        })
        .collect()
}

fn global_assess(answer: &HashMap<u64, usize>, dataset: &DataSet) -> std::io::Result<()> {
    let mut counts: Vec<HashMap<_, usize>> = (0..3).map(|_| HashMap::new()).collect();
    for asn in dataset.assignments.iter() {
        let ans = answer[&asn.id];
        *counts[ans].entry(asn.cluster).or_default() += 1;
    }
    let mut out = BufWriter::new(std::fs::File::create("global.tsv")?);
    writeln!(&mut out, "Answer\tPred\tCount")?;
    for (ans, preds) in counts.iter().enumerate() {
        let mut preds: Vec<_> = preds.into_iter().collect();
        preds.sort_by_key(|x| x.0);
        for (pred, count) in preds {
            let ans = match ans {
                0 => "hapA",
                1 => "hapB",
                2 => "hapC",
                _ => unreachable!(),
            };
            writeln!(&mut out, "{}\t{}\t{}", ans, pred, count)?;
        }
    }
    Ok(())
}

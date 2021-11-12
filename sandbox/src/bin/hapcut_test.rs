use definitions::*;
use haplotyper::Encode;
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    use haplotyper::hapcut::HapCut;
    use haplotyper::hapcut::HapCutConfig;
    let copy_number: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|x| (x.id, x.cluster_num))
        .collect();
    let seed = 342;
    let thr = 0.1;
    let flip = 0.8;
    // for seed in vec![342094, 4234, 56, 43290] {
    //     for thr in vec![0.05, 0.1, 0.15, 0.2, 0.3] {
    let config = HapCutConfig::new(seed, thr, flip, 0, copy_number.clone());
    // println!("{}\t{}", seed, thr);
    let (_consis, phases) = ds.hapcut(&config);
    let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    {
        let phases: HashMap<(u64, u64), _> = phases.iter().flat_map(|p| p.phases()).collect();
        let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
        for read in ds.encoded_reads.iter() {
            let ans = id2name[&read.id].contains("000251v2") as usize;
            for node in read.nodes.iter() {
                counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
            }
        }
        let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
        println!("UNIT\tunit\tcluster\thap1\thap2\tphase\tpurity\tscore");
        for (&(unit, cluster), counts) in counts.iter() {
            let phase = phases[&(unit, cluster)];
            let score = score[&unit];
            let total = counts[0] + counts[1];
            let pur = counts[0].max(counts[1]) as f64 / total as f64;
            println!(
                "UNIT\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
                unit, cluster, counts[0], counts[1], phase, pur, score
            );
        }
    }
    {
        let mut result: HashMap<_, [u32; 2]> = HashMap::new();
        for &Assignment { id, cluster } in ds.assignments.iter() {
            *result
                .entry(cluster)
                .or_default()
                .get_mut(id2name[&id].contains("251v2") as usize)
                .unwrap() += 1;
        }
        println!("PB\tPhaseBlockID\tReads from DBB\tReads from QBL");
        let mut result: Vec<_> = result.iter().collect();
        result.sort_by_key(|x| x.0);
        for (cluster, haps) in result.iter() {
            println!("PB\t{}\t{}\t{}", cluster, haps[0], haps[1]);
        }
    }
    // for &Assignment { id, cluster } in ds.assignments.iter() {
    //     *result
    //         .entry(cluster)
    //         .or_default()
    //         .get_mut(answer[&id] as usize)
    //         .unwrap() += 1;
    // }

    // eprintln!("{}", ds.raw_reads.len());
    // let total: HashSet<_> = ds.assignments.iter().map(|x| x.id).collect();
    // let unphased: usize = ds
    //     .raw_reads
    //     .iter()
    //     .filter(|r| !total.contains(&r.id))
    //     .count();
    // eprintln!("{}\t{}", unphased, total.len());
    // let max_cl = ds.assignments.iter().map(|x| x.cluster).max().unwrap();
    // let assignments: HashMap<_, usize> = ds.assignments.iter().map(|a| (a.id, a.cluster)).collect();
    // for cl in 0..=max_cl {
    //     let mut counts: HashMap<_, HashMap<_, u32>> = HashMap::new();
    //     for read in ds
    //         .encoded_reads
    //         .iter()
    //         .filter(|r| assignments.get(&r.id) == Some(&cl))
    //     {
    //         for node in read.nodes.iter() {
    //             *counts
    //                 .entry(node.unit)
    //                 .or_default()
    //                 .entry(node.cluster)
    //                 .or_default() += 1;
    //         }
    //     }
    //     for (unit, counts) in counts.iter() {
    //         let mut counts: Vec<_> = counts.iter().map(|(c, x)| (c, x)).collect();
    //         counts.sort_by_key(|x| x.0);
    //         for (cluster, count) in counts.iter() {
    //             println!("{}\t{}\t{}\t{}", cl, unit, cluster, count);
    //         }
    //     }
    // }
    //     }
    // }
    /////

    Ok(())
}

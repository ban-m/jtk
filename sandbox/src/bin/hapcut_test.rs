use definitions::*;
use haplotyper::Encode;
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
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
    let seed = 34220392;
    let thr = 0.8;
    let flip = 0.8;
    // for seed in vec![342094, 4234, 56, 43290] {
    //     for thr in vec![0.05, 0.1, 0.15, 0.2, 0.3] {
    let config = HapCutConfig::new(seed, thr, flip, 0, 2, copy_number.clone());
    // println!("{}\t{}", seed, thr);
    // Want to check the correct phasing -> OK
    // Collapse errorneous variants -> NG.
    // let correct_phases = {
    //     ds = squish_errorneous_clustering(ds);
    //     // Create the correct phaseing.
    //     haplotyper::local_clustering::normalize::normalize_local_clustering(&mut ds);
    // let phases: Vec<_> = predict_correct_phasing(&ds);
    //     let reads: Vec<_> = ds
    //         .encoded_reads
    //         .iter()
    //         .map(|read| {
    //             let nodes: Vec<_> = read.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
    //             (read.id, nodes)
    //         })
    //         .collect();
    //     let consis = haplotyper::hapcut::Phase::consistency(&phases, &reads, &config);
    //     eprintln!("ConsistencyOfCorrectPhase\t{}", consis);
    //     phases
    // };
    let sq = ds.hapcut_squish_diplotig(&config);
    eprintln!("HAPCUT\tSquish\t{:?}", sq);
    env_logger::init();
    let (_consis, phases) = ds.hapcut(&config);
    // {
    //     let mut counts: HashMap<_, u32> = HashMap::new();
    //     for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
    //         *counts.entry(node.unit).or_default() += 1;
    //     }
    //     assert_eq!(correct_phases.len(), phases.len());
    //     for (unit, (correct, pred)) in correct_phases.iter().zip(phases.iter()).enumerate() {
    //         if correct != pred {
    //             eprintln!(
    //                 "{}\t{}\t{:?}\t{:?}",
    //                 unit,
    //                 counts[&(unit as u64)],
    //                 correct.phase(),
    //                 pred.phase(),
    //             );
    //         }
    //     }
    // }
    let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    {
        let phases: HashMap<(u64, u64), _> = phases.iter().flat_map(|p| p.phases()).collect();
        let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
        for read in ds.encoded_reads.iter() {
            //let ans = id2name[&read.id].contains("000251v2") as usize;
            let ans = id2name[&read.id].contains("hapA") as usize;
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

fn squish_errorneous_clustering(mut ds: DataSet) -> DataSet {
    let max_cluster_num = get_max_cluster_num(&ds);
    let counts = get_counts(&ds);
    let bad_clusterings: HashSet<_> = counts
        .iter()
        .filter(|(_, counts)| {
            let frac = counts[0].max(counts[1]) as f64 / (counts[0] + counts[1]) as f64;
            frac < 0.8
        })
        .map(|(&x, _)| x)
        .collect();
    let bad_unit: HashSet<u64> = max_cluster_num
        .iter()
        .filter(|(&unit, &clsize)| (0..clsize).all(|cl| bad_clusterings.contains(&(unit, cl))))
        .map(|(&unit, _)| unit)
        .collect();
    eprintln!("Squish {} units.", bad_unit.len());
    for read in ds.encoded_reads.iter_mut() {
        for node in read.nodes.iter_mut() {
            if bad_unit.contains(&node.unit) {
                node.cluster = 0;
            }
            if node.unit == 1625 {
                node.cluster = 0;
            }
        }
    }
    ds
}

fn predict_correct_phasing(ds: &DataSet) -> Vec<haplotyper::hapcut::Phase> {
    let max_cluster_num = get_max_cluster_num(ds);
    let counts = get_counts(ds);
    let max_num_unit: u64 = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    (0..=max_num_unit)
        .map(|unit| {
            let cl_num = max_cluster_num.get(&unit).unwrap_or(&0);
            let phase: Vec<i8> = (0..*cl_num)
                .map(|cl| {
                    counts
                        .get(&(unit, cl))
                        .map(|counts| {
                            let total = counts[0] + counts[1];
                            if total * 7 / 10 < counts[0] {
                                1
                            } else if total * 7 / 10 < counts[1] {
                                -1
                            } else {
                                0
                            }
                        })
                        .unwrap_or(0)
                })
                .collect();
            haplotyper::hapcut::Phase::new(unit, phase)
        })
        .collect()
}

fn get_max_cluster_num(ds: &DataSet) -> HashMap<u64, u64> {
    let mut max_cluster_num: HashMap<_, u64> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        max_cluster_num
            .entry(node.unit)
            .and_modify(|x| *x = (*x).max(node.cluster + 1))
            .or_insert(node.cluster + 1);
    }
    max_cluster_num
}

fn get_counts(ds: &DataSet) -> HashMap<(u64, u64), [u32; 2]> {
    let id2name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = id2name[&read.id].contains("000251v2") as usize;
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    counts
}

// use kiley::bialignment::edit_dist;
// use rand::Rng;
use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use log::*;
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.name.clone()))
        .collect();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(3214);
    let unit = 429;
    ////////////////
    // let seqs: Vec<_> = ds
    //     .encoded_reads
    //     .iter()
    //     .flat_map(|read| {
    //         let is_hap1 = id2desc[&read.id].starts_with("hapA");
    //         read.nodes
    //             .iter()
    //             .filter(|n| n.unit == unit)
    //             .map(|x| (is_hap1, x.seq()))
    //             .collect::<Vec<_>>()
    //     })
    //     .collect();
    // let all_seq: Vec<_> = seqs.iter().map(|x| x.1).collect();
    // use rand::Rng;
    // debug!("Ternary");
    // let consed = kiley::ternary_consensus(&all_seq, rng.gen(), 3, 50);
    // use kiley::gphmm::*;
    // let model = GPHMM::<Cond>::new_three_state(0.9, 0.05, 0.10, 0.95);
    // let consed = model
    //     .correct_until_convergence_banded(&consed, &all_seq, 50)
    //     .unwrap();
    // let mut hap_a: Vec<_> = seqs.iter().filter_map(|x| x.0.then(|| x.1)).collect();
    // hap_a.extend(all_seq.iter().take(30 - hap_a.len().min(30)));
    // debug!("Polishing A.");
    // let cons_a = model
    //     .correct_until_convergence_banded(&consed, &hap_a, 50)
    //     .unwrap();
    // let mut hap_c: Vec<_> = seqs.iter().filter_map(|x| (!x.0).then(|| x.1)).collect();
    // hap_c.extend(all_seq.iter().take(30 - hap_c.len().min(30)));
    // debug!("Polishing C.");
    // let cons_c = model
    //     .correct_until_convergence_banded(&consed, &hap_c, 50)
    //     .unwrap();
    // debug!("Polished");
    // let (dist, op) = kiley::bialignment::edit_dist_slow_ops(&cons_a, &cons_c);
    // println!("{}", dist);
    // let (ar, opr, cr) = kiley::bialignment::recover(&cons_a, &cons_c, &op);
    // for ((ar, opr), cr) in ar.chunks(200).zip(opr.chunks(200)).zip(cr.chunks(200)) {
    //     println!("{}", String::from_utf8_lossy(ar));
    //     println!("{}", String::from_utf8_lossy(opr));
    //     println!("{}\n", String::from_utf8_lossy(cr));
    // }
    // for &(is_a, seq) in seqs.iter() {
    //     let hapa = edlib_sys::global_dist(seq, &cons_a);
    //     let hapc = edlib_sys::global_dist(seq, &cons_c);
    //     println!("SCORE\t{}\t{}", is_a as u8, hapa as i32 - hapc as i32);
    //     println!("SEQ\t{}\t{}\n", is_a as u8, String::from_utf8_lossy(seq));
    // }
    //////////////////
    let chunks: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|read| read.nodes.iter().filter(|n| n.unit == unit))
        .cloned()
        .collect();
    let seqs: Vec<_> = chunks.iter().map(|x| x.seq()).collect();
    let answers: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|read| {
            let answer = id2desc[&read.id].starts_with("hapA");
            let len = read.nodes.iter().filter(|n| n.unit == unit).count();
            vec![answer as u8; len]
        })
        .collect();
    let config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, 20, 5, 2, 20);
    debug!("Start");
    let asn = haplotyper::local_clustering::kmeans::clustering(&seqs, &mut rng, &config);
    for (&asn, answer) in asn.iter().zip(answers.iter()) {
        println!("DUMP\t{}\t{}\tKMEANS", asn, answer,);
    }
    ////////////////
    let mut config = haplotyper::ClusteringConfig::default();
    config.cluster_num = 2;
    let mut ds = ds;
    let mut chunks: Vec<_> = ds
        .encoded_reads
        .iter_mut()
        .flat_map(|read| {
            let id = read.id;
            read.nodes
                .iter_mut()
                .filter(|n| n.unit == unit)
                .map(|n| (id, 0, n))
                .collect::<Vec<_>>()
        })
        .collect();
    let refunit = ds.selected_chunks.iter().find(|u| u.id == unit).unwrap();
    haplotyper::unit_clustering(&mut chunks, &config, refunit, 1);
    for (i, ((id, _, node), asn)) in chunks.iter().zip(asn.iter()).enumerate() {
        println!(
            "DUMP\t{}\t{}\t{}\t{}\t{}",
            i, id, node.cluster, asn, id2desc[id]
        );
    }
    let answer: Vec<_> = chunks
        .iter()
        .map(|(id, _, _)| id2desc[id].starts_with("hapA") as u8)
        .collect();
    let km = haplotyper::local_clustering::rand_index(&answer, &asn);
    let pred: Vec<_> = chunks.iter().map(|(_, _, n)| n.cluster as u8).collect();
    let default = haplotyper::local_clustering::rand_index(&answer, &pred);
    eprintln!("{},{}", km, default);
    Ok(())
}

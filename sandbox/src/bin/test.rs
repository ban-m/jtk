// use definitions::*;
// use poa_hmm::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
// use std::io::BufReader;

fn main() {
    let len = 2000;
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(342389);
    let template = kiley::gen_seq::generate_seq(&mut rng, len);
    let hap_div = kiley::gen_seq::Profile {
        sub: 0.001 / 3f64,
        ins: 0.001 / 3f64,
        del: 0.001 / 3f64,
    };
    let haps = vec![
        kiley::gen_seq::introduce_randomness(&template, &mut rng, &hap_div),
        kiley::gen_seq::introduce_randomness(&template, &mut rng, &hap_div),
    ];
    println!(
        "{}",
        kiley::alignment::bialignment::edit_dist(&haps[0], &haps[1])
    );
    let coverage = 20;
    let num_hap = 2;
    let prof = kiley::gen_seq::PROFILE;
    let xs: Vec<_> = (0..coverage * num_hap)
        .map(|i| {
            let index = i / coverage;
            kiley::gen_seq::introduce_randomness(&haps[index], &mut rng, &prof)
        })
        .collect();
    let mut assignments: Vec<_> = xs.iter().map(|_| rng.gen_range(0..num_hap)).collect();
    for t in 0..10 {
        let consensus: Vec<_> = (0..num_hap)
            .map(|hap| {
                let xs: Vec<_> = xs
                    .iter()
                    .zip(assignments.iter())
                    .filter_map(|(x, &asn)| (asn == hap).then(|| x.as_slice()))
                    .collect();
                kiley::consensus_poa(&xs, rng.gen(), 10, 10, "CLR")
            })
            .collect();
        for (asn, x) in assignments.iter_mut().zip(xs.iter()) {
            *asn = consensus
                .iter()
                .enumerate()
                .map(|(idx, cns)| (idx, kiley::alignment::bialignment::edit_dist(x, &cns)))
                .min_by_key(|x| x.1)
                .map(|x| x.0)
                .unwrap();
        }
        println!("{}-th iter ended.", t);
        use std::collections::HashMap;
        let mut counts: HashMap<_, u32> = HashMap::new();
        for &asn in assignments.iter() {
            *counts.entry(asn).or_default() += 1;
        }
        println!("{:?}", counts);
    }
    for (idx, asn) in assignments.iter().enumerate() {
        println!("{},{}", idx, asn);
    }
    // let args: Vec<_> = std::env::args().collect();
    // let ds: DataSet =
    //     serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
    //         .unwrap();
    // let id2name: std::collections::HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|read| (read.id, read.name.to_string()))
    //     .collect();
    // for read in ds.encoded_reads.iter() {
    //     for edge in read.edges.iter() {
    //         if edge.from == edge.to {
    //             println!("{}\t{:?}", id2name[&read.id], edge);
    //         }
    //     }
    // }
}

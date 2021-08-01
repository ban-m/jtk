use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use poa_hmm::gen_sample;
    use rand::Rng;
    let seed = 923480;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let len = 10;
    let test_num = 40;
    let former_half = gen_sample::generate_seq(&mut rng, len);
    let latter_half = gen_sample::generate_seq(&mut rng, len);
    let rep = b"AT";
    // let rep = b"A";
    let clusters = 2;
    let mut templates = vec![];
    for i in 0..clusters {
        let rep_num = 5 + 3 * i;
        // let rep_range = 10..20;
        // let rep_num = rng.gen_range(rep_range.clone());
        log::debug!("{}\t{}", i, rep_num);
        let mut template = former_half.clone();
        for _ in 0..rep_num {
            template.extend_from_slice(rep);
        }
        template.extend_from_slice(&latter_half);
        templates.push(template);
    }
    let profile = gen_sample::PROFILE.norm().mul(0.15);
    let probs = vec![0.5, 0.5];
    use sandbox::generate_test_data;
    let (dataset, answer) = generate_test_data(&templates, test_num, &mut rng, &probs, &profile);
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let all_patterns: Vec<_> = (5..9)
        .map(|rep_num| {
            let mut buf = former_half.clone();
            for _ in 0..rep_num {
                buf.extend_from_slice(rep);
            }
            buf.extend_from_slice(&latter_half);
            buf
        })
        .collect();
    for (read, asn) in dataset.iter().zip(answer) {
        let all_patterns: Vec<_> = all_patterns
            .iter()
            .map(|x| hmm.likelihood(x, read))
            .collect();
        let all_patterns: Vec<_> = all_patterns
            .iter()
            .map(|x| x - all_patterns[0])
            .map(|x| format!("{:.2}", x))
            .collect();
        let read = String::from_utf8_lossy(read);
        println!("{}\t{}\t{}", asn, all_patterns.join("\t"), read);
    }
    // let args: Vec<_> = std::env::args().collect();
    // let start = std::time::Instant::now();
    // let mut ds: DataSet =
    //     serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    // let end = std::time::Instant::now();
    // eprintln!("Open\t{}", (end - start).as_secs());
    // use std::collections::HashMap;
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // use haplotyper::em_correction::ClusteringCorrection;
    // let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    // let ds = ds.correct_clustering_em(5, 5, 3);

    // let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     // let ans = id2desc[&read.id].contains("hapA") as usize;
    //     let ans = id2desc[&read.id].contains("000252v2") as usize;
    //     for node in read.nodes.iter() {
    //         counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
    //     }
    // }
    // for ((unit, cluster), counts) in counts.iter() {
    //     println!("{}\t{}\t{}\t{}", unit, cluster, counts[0], counts[1]);
    // }
    Ok(())
}

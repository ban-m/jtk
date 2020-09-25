use std::collections::HashMap;
fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let sam: Vec<_> = bio_utils::sam::load_sam_file(&args[1]).unwrap();
    let reference: HashMap<_, _> = bio_utils::fasta::parse_into_vec(&args[2])
        .unwrap()
        .into_iter()
        .map(|r| (r.id().to_string(), r))
        .collect();
    println!("RefName\tMat\tMism\tIns\tDel");
    for record in sam {
        let (mut ins, mut del, mut mat, mut mism) = (0, 0, 0, 0);
        use bio_utils::sam::Op;
        for op in record.cigar() {
            match op {
                Op::Insertion(l) => ins += l,
                Op::Deletion(l) => del += l,
                Op::Match(l) => mat += l,
                Op::Align(l) => mat += l,
                Op::Mismatch(l) => mism += l,
                _ => {}
            }
        }
        // let (mut max_ins, mut max_del) = (0, 0);
        // for op in record.cigar() {
        //     match op {
        //         Op::Insertion(l) => max_ins = max_ins.max(l),
        //         Op::Deletion(l) => max_del = max_del.max(l),
        //         _ => {}
        //     }
        // }
        // println!("{}\t{}", max_ins, max_del);
        println!("{}\t{}\t{}\t{}\t{}", record.r_name(), mat, mism, ins, del);
        let seq = record.seq().as_bytes();
        if seq.len() > 1 {
            let cigar = record.cigar();
            let refr = reference[record.r_name()].seq();
            let (q, o, r) = bio_utils::sam::recover_alignment(&cigar, seq, &refr, record.pos());
            for ((q, o), r) in q.chunks(150).zip(o.chunks(150)).zip(r.chunks(150)) {
                println!("{}", String::from_utf8_lossy(q));
                println!("{}", String::from_utf8_lossy(o));
                println!("{}\n", String::from_utf8_lossy(r));
            }
        }
    }
    // for i in 0..100 {
    //     let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(i as u64);
    //     let template = poa_hmm::gen_sample::generate_seq(&mut rng, 100);
    //     let p = &poa_hmm::gen_sample::PROFILE;
    //     let seqs: Vec<_> = (0..10)
    //         .map(|_| poa_hmm::gen_sample::introduce_randomness(&template, &mut rng, p))
    //         .collect();
    //     let seqs: Vec<_> = (0..10)
    //         .map(|_| {
    //             let mut seqs: Vec<_> = seqs.clone();
    //             use rand::seq::SliceRandom;
    //             seqs.shuffle(&mut rng);
    //             poa_hmm::POA::from_vec_default(&seqs).consensus_homopolymer()
    //         })
    //         .collect();
    //     let poa = poa_hmm::POA::from_vec_default(&seqs);
    //     let pred = poa.consensus_homopolymer();
    //     let dist1 = bio_utils::alignments::edit_dist(&template, &pred);
    //     // println!("{}", String::from_utf8_lossy(&template));
    //     // println!("{}", String::from_utf8_lossy(&pred));
    //     println!("{}\t{}", i, dist1);
    // }
}

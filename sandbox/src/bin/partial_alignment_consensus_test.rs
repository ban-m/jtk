use kiley::gen_seq;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let len = 2000;
    let seed = 2032;
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    let seq = gen_seq::generate_seq(&mut rng, len);
    let coverage = 20;
    let profile = gen_seq::Profile {
        sub: 0.02,
        del: 0.02,
        ins: 0.02,
    };
    let draft = gen_seq::introduce_randomness(&seq, &mut rng, &profile);
    let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let length = rng.gen_range(4 * len / 5..len);
    let template = match i % 2 == 0 {
        true => &seq[..length],
        false => &seq[len - length..],
    };
    let query = gen_seq::introduce_randomness(template, &mut rng, &profile);
    hmm.modification_table(&draft, &query, 100, &ops);
    // let full_lengths: Vec<_> = (0..coverage)
    //     .map(|_| gen_seq::introduce_randomness(&seq, &mut rng, &profile))
    //     .collect();
    // let partial: Vec<_> = (0..coverage)
    //     .map(|i| {
    //         let length = rng.gen_range(4 * len / 5..len);
    //         let template = match i % 2 == 0 {
    //             true => &seq[..length],
    //             false => &seq[len - length..],
    //         };
    //         gen_seq::introduce_randomness(template, &mut rng, &profile)
    //     })
    //     .collect();
    // use kiley::op::Op;
    // let op2op = [Op::Match, Op::Ins, Op::Del, Op::Mismatch];
    // let (seqs, mut ops): (Vec<_>, Vec<_>) = full_lengths
    //     .iter()
    //     .chain(partial.iter())
    //     .map(|query| {
    //         let ops = edlib_sys::global(&draft, query);
    //         let ops: Vec<_> = ops.into_iter().map(|op| op2op[op as usize]).collect();
    //         hmm.modification_table(&draft, &query, 100, &ops);
    //         (query.as_slice(), ops)
    //     })
    //     .unzip();
    // println!("kfjds");
    // let consensus = hmm.polish_until_converge_with(&draft, &seqs, &mut ops, 100);
    // let dist = edlib_sys::global_dist(&seq, &consensus);
    // let seq = std::str::from_utf8(&seq).unwrap();
    // let consensus = std::str::from_utf8(&consensus).unwrap();
    // println!("{seq}\n{consensus}\n{dist}");
    Ok(())
}

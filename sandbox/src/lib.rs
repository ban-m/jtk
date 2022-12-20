pub const IS_MOCK: bool = false;
use rand::Rng;
pub fn test() {
    let xs: Vec<_> = (0..1000).collect();
    let i = 32948;
    for x in xs[i..].iter() {
        println!("{}", x);
    }
    for x in xs[..i].iter() {
        println!("{}", x);
    }
}

pub fn generate_test_data<T: Rng>(
    templates: &[Vec<u8>],
    test_num: usize,
    rng: &mut T,
    probs: &[f64],
    profile: &kiley::gen_seq::Profile,
) -> (Vec<Vec<u8>>, Vec<usize>) {
    let choices: Vec<_> = (0..templates.len()).collect();
    use rand::seq::SliceRandom;
    let mut answer: Vec<_> = (0..test_num)
        .filter_map(|_| choices.choose_weighted(rng, |&k| probs[k]).ok().copied())
        .collect();
    answer.sort_unstable();
    let dataset: Vec<_> = answer
        .iter()
        .map(|&idx| kiley::gen_seq::introduce_randomness(&templates[idx], rng, profile))
        .collect();
    assert_eq!(dataset.len(), answer.len());
    (dataset, answer)
}

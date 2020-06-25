use rand_distr::{Distribution, Normal};
fn main() {
    let normal = Normal::new(100., 10.).unwrap();
    let mut rng = rand::thread_rng();
    let data: Vec<_> = normal
        .sample_iter(&mut rng)
        .take(10000)
        .map(|x: f64| if x < 0. { 0 } else { x.floor() as usize })
        .collect();
    let hst = histgram_viz::Histgram::new(&data);
    println!("{}", hst.format(40, 10));
}

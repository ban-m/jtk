type FeatureVector = (
    Vec<Vec<f64>>,
    Vec<(usize, crate::likelihood_gains::DiffType)>,
);
type ClusteringDevResult = (Vec<usize>, Vec<Vec<f64>>, f64, usize);

pub fn cluster_filtered_variants_exact(
    (variants, variant_type): &FeatureVector,
    copy_num: usize,
) -> ClusteringDevResult {
    let feature_dim = variant_type.len();
    let mut selected_variants: Vec<_> = vec![0; copy_num];
    let choises = 1 << feature_dim;
    let last_loop = vec![choises - 1; copy_num];
    let mut max = 0f64;
    let mut argmax = get_result(&selected_variants, variants);
    while selected_variants != last_loop {
        let score = calc_score(&selected_variants, variants);
        if max < score {
            argmax = get_result(&selected_variants, variants);
            max = score;
        }
        increment_one(&mut selected_variants, choises);
    }
    argmax
}

fn get_result(vars: &[usize], variants: &[Vec<f64>]) -> ClusteringDevResult {
    let score = calc_score(vars, variants);
    let (assignments, lk_gain): (Vec<_>, Vec<_>) = variants
        .iter()
        .map(|xs| {
            let lk_gain: Vec<_> = vars
                .iter()
                .map(|&selection| get_exact_score(selection, xs))
                .collect();
            let (max_id, _) = lk_gain
                .iter()
                .enumerate()
                .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
                .unwrap();
            (max_id, lk_gain)
        })
        .unzip();
    (assignments, lk_gain, score, vars.len())
}

fn get_exact_score(selection: usize, xs: &[f64]) -> f64 {
    xs.iter()
        .enumerate()
        .filter_map(|(i, x)| (((1 << i) & selection) != 0).then_some(x))
        .sum()
}

fn calc_score(vars: &[usize], variants: &[Vec<f64>]) -> f64 {
    fn max_gain(vars: &[usize], xs: &[f64]) -> f64 {
        vars.iter()
            .map(|&selection| get_exact_score(selection, xs))
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }
    variants.iter().map(|xs| max_gain(vars, xs)).sum()
}

fn increment_one(vars: &mut [usize], max: usize) {
    let mut idx = 0;
    while max == vars[idx] + 1 {
        idx += 1;
    }
    vars[idx] += 1;
    for j in 0..idx {
        vars[j] = vars[idx];
    }
    for w in vars.windows(2) {
        assert!(w[1] <= w[0]);
    }
}

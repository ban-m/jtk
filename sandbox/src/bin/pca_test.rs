use bio::alignment::pairwise::Aligner;
use log::*;
use nalgebra::*;
use poa_hmm::gen_sample;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    debug!("Start");
    let cluster_num = 2;
    let len = 50;
    let chain = 1;
    let seed = 923;
    let num = 20;
    let _div = gen_sample::Profile {
        sub: 0.01,
        ins: 0.01,
        del: 0.01,
    };
    let p = gen_sample::PROFILE;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let template: Vec<_> = (0..chain)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates: Vec<Vec<_>> = (0..cluster_num)
        .map(|_| {
            template
                .iter()
                //.map(|t| gen_sample::introduce_randomness(t, &mut rng, &div))
                .map(|t| gen_sample::introduce_errors(t, &mut rng, 1, 0, 0))
                .collect()
        })
        .collect();
    for (i, t) in templates.iter().enumerate() {
        for (j, s) in templates.iter().enumerate().skip(i + 1) {
            let dist = t
                .iter()
                .zip(s)
                .map(|(a, b)| bio::alignment::distance::levenshtein(a, b))
                .sum::<u32>();
            eprintln!("{}\t{}\t{}", i, j, dist);
        }
    }
    let reads = gen_reads(&templates, num, &mut rng, &p);
    debug!("Gen reads:{}", reads.len());
    let consensus = take_consensus(&reads, chain);
    debug!("Gen consensus:{}", consensus.len());
    let alns = convert_to_alignment(&reads, &consensus);
    let summary = get_summary(&alns, alns.iter().map(|(_, xs)| xs.len()).max().unwrap());
    for (i, _) in summary.iter().enumerate().filter(|x| x.1.is_some()) {
        println!("{}", i);
    }
    for (id, aln) in alns.iter() {
        let mut x: String = String::new();
        for c in aln.iter() {
            let c = match *c {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                4 => '-',
                _ => continue,
            };
            x.push(c);
        }
        println!("{:3} {}", id, x);
    }
    {
        let mut x = String::new();
        for cs in consensus.iter() {
            for &c in cs.iter() {
                x.push(c as char);
                x.push('-');
            }
        }
        println!("{:3} {}", 0, x)
    };
    debug!("Aligned reads");
    let features = principal_component_analysis(&alns, 4);
    debug!("Done PCA");
    for (idx, features) in features {
        let features: Vec<_> = features.iter().map(|x| format!("{}", x)).collect();
        println!("{},{}", idx, features.join(","));
    }
}

fn get_summary(alns: &[(usize, Vec<u8>)], len: usize) -> Vec<Option<(u8, u8)>> {
    let mut count = vec![vec![0; 5]; len];
    for (_, s) in alns.iter() {
        for (i, &x) in s.iter().enumerate() {
            count[i][x as usize] += 1;
        }
    }
    count
        .into_iter()
        .map(|cs| {
            let (argmax, max) = cs.iter().enumerate().max_by_key(|x| x.1).unwrap();
            let (next, next_max) = cs
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != argmax)
                .max_by_key(|x| x.1)
                .unwrap();
            if *max < 3 * *next_max {
                Some((argmax as u8, next as u8))
            } else {
                None
            }
        })
        .collect()
}

fn principal_component_analysis(alns: &[(usize, Vec<u8>)], dim: usize) -> Vec<(usize, Vec<f64>)> {
    let len = alns.iter().map(|(_, xs)| xs.len()).max().unwrap();
    let summary = get_summary(alns, len);
    let len = summary.iter().filter(|x| x.is_some()).count();
    debug!("Length:{}", len);
    let vectors: Vec<_> = alns
        .iter()
        .map(|(_, xs)| {
            let xs = xs.iter().zip(summary.iter()).filter_map(|(&x, info)| {
                info.map(|(m, n)| {
                    if x == m {
                        1.
                    } else if x == n {
                        -1.
                    } else {
                        0.
                    }
                })
            });
            DVector::from_iterator(len, xs)
        })
        .collect();
    for (idx, v) in vectors.iter().enumerate() {
        let v: Vec<_> = v.iter().map(|v| format!("{:.0}", v)).collect();
        eprintln!("{}\t{}", idx, v.join("\t"));
    }
    let mean = vectors.iter().fold(DVector::zeros(len), |x, y| x + y) / vectors.len() as f64;
    debug!("Calc mean");
    let covariance = vectors
        .iter()
        .map(|x| (x - mean.clone()) * (x - mean.clone()).transpose())
        .fold(DMatrix::zeros(len, len), |x, y| x + y)
        / vectors.len() as f64;
    debug!("Calc coval");
    let eigens = covariance.clone().symmetric_eigen();
    let mut eigen_and_eigenvec: Vec<_> = eigens
        .eigenvectors
        .column_iter()
        .zip(eigens.eigenvalues.iter())
        .collect();
    eigen_and_eigenvec.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
    eigen_and_eigenvec.reverse();
    debug!("Calc Eigen Values");
    let pca_vectors = &eigen_and_eigenvec[..dim];
    vectors
        .iter()
        .zip(alns)
        .map(|(x, &(i, _))| {
            let features: Vec<_> = pca_vectors.iter().map(|(v, _)| x.dot(&v)).collect();
            (i, features)
        })
        .collect()
}

fn gen_reads<R: rand::Rng>(
    templates: &[Vec<Vec<u8>>],
    num: usize,
    rng: &mut R,
    p: &gen_sample::Profile,
) -> Vec<(usize, Vec<Vec<u8>>)> {
    templates
        .iter()
        .enumerate()
        .flat_map(|(idx, t)| {
            (0..num)
                .map(|_| {
                    let seq: Vec<_> = t
                        .iter()
                        .map(|s| gen_sample::introduce_randomness(s, rng, p))
                        .collect();
                    (idx, seq)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn take_consensus(reads: &[(usize, Vec<Vec<u8>>)], chain: usize) -> Vec<Vec<u8>> {
    let mut chunks = vec![vec![]; chain];
    for (_, r) in reads.iter() {
        for (i, c) in r.iter().enumerate() {
            chunks[i].push(c.as_slice());
        }
    }
    chunks
        .iter()
        .map(|cs| poa_hmm::POA::from_slice_default(cs).consensus())
        .collect()
}

fn convert_to_alignment(
    reads: &[(usize, Vec<Vec<u8>>)],
    consensus: &[Vec<u8>],
) -> Vec<(usize, Vec<u8>)> {
    let mut aligner = Aligner::new(-3, -1, |x, y| if x == y { 2 } else { -5 });
    reads
        .iter()
        .map(|(idx, r)| {
            let aln: Vec<_> = r
                .iter()
                .enumerate()
                .flat_map(|(idx, seq)| alignment(&mut aligner, seq, &consensus[idx]))
                .collect();
            (*idx, aln)
        })
        .collect()
}

fn alignment<F: bio::alignment::pairwise::MatchFunc>(
    aligner: &mut Aligner<F>,
    query: &[u8],
    template: &[u8],
) -> Vec<u8> {
    let aln = aligner.global(&query, &template);
    let mut seq = vec![];
    let mut ins = vec![b'-'; template.len()];
    let (mut rpos, mut qpos) = (0, 0);
    let mut prev = None;
    for &op in aln.operations.iter() {
        use bio::alignment::AlignmentOperation::*;
        match op {
            Del => {
                seq.push(b'-');
                rpos += 1;
            }
            Ins => {
                if prev != Some(Ins) && rpos < ins.len() {
                    ins[rpos] = query[qpos];
                }
                qpos += 1;
            }
            Subst | Match => {
                seq.push(query[qpos]);
                rpos += 1;
                qpos += 1;
            }
            _ => panic!(),
        }
        prev = Some(op);
    }
    let mut result = vec![];
    for (x, y) in seq.into_iter().zip(ins) {
        result.push(x);
        result.push(y);
    }
    result
        .iter()
        .map(|&b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'-' => 4,
            _ => panic!(),
        })
        .collect()
}

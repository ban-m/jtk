// use bio::alignment::pairwise::Aligner;
// use nalgebra::*;
// use rand::Rng;
// use rand_xoshiro::Xoshiro256PlusPlus;
use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let dataset: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new).unwrap())
            .unwrap();
    for read in dataset.encoded_reads.iter() {
        let mut units: HashMap<_, u32> = HashMap::new();
        for node in read.nodes.iter() {
            *units.entry(node.unit).or_default() += 1;
        }
        if units.values().any(|&x| x > 1) {
            let r: Vec<_> = read.nodes.iter().map(|n| format!("{}", n.unit)).collect();
            eprintln!("{}", r.join("-"));
        }
        for (&unit_id, _) in units.iter().filter(|&(_, &val)| val > 1) {
            let unit = dataset
                .selected_chunks
                .iter()
                .find(|u| u.id == unit_id)
                .unwrap();
            for node in read.nodes.iter() {
                eprintln!("DUMP Alignment {}", unit_id);
                let (q, a, r) = node.recover(unit);
                for ((q, a), r) in q.chunks(100).zip(a.chunks(100)).zip(r.chunks(100)) {
                    eprintln!("{}", String::from_utf8_lossy(q));
                    eprintln!("{}", String::from_utf8_lossy(a));
                    eprintln!("{}", String::from_utf8_lossy(r));
                    eprintln!();
                }
            }
        }
    }
    // let num_encoded = dataset.encoded_reads.len();
    // let num_nodes = dataset
    //     .encoded_reads
    //     .iter()
    //     .map(|er| er.nodes.len())
    //     .sum::<usize>();
    // eprintln!("{}", num_nodes);
    // eprintln!(
    //     "{} encodes out of {}, {} units.",
    //     dataset.raw_reads.len(),
    //     num_encoded,
    //     num_nodes
    // );
    // let mut counts: HashMap<_, usize> = HashMap::new();
    // for read in dataset.encoded_reads.iter() {
    //     for node in read.nodes.iter() {
    //         *counts.entry(node.unit).or_default() += 1;
    //     }
    // }
    // let mut counts: Vec<_> = counts.into_iter().collect();
    // counts.sort_by_key(|x| x.0);
    // for (u, v) in counts.iter() {
    //     eprintln!("{}\t{}", u, v);
    // }
    Ok(())
}

//     let args: Vec<_> = std::env::args().collect();
//     let mut buf = String::new();
//     let mut rdr = BufReader::new(File::open(&args[1]).unwrap());
//     rdr.read_to_string(&mut buf).unwrap();
//     let reads: Vec<_> = buf
//         .split('>')
//         .filter_map(|record| {
//             let mut record = record.split('\n');
//             let name = record.next().unwrap();
//             let seq: Vec<_> = record.flat_map(|x| x.as_bytes()).copied().collect();
//             if seq.is_empty() {
//                 None
//             } else {
//                 Some((name, seq))
//             }
//         })
//         .collect();
//     let queries: Vec<_> = reads.iter().map(|x| x.1.as_slice()).collect();
//     let template = poa_hmm::POA::from_slice_default(&queries).consensus();
//     let mut aligner = Aligner::new(-3, -1, |x, y| if x == y { 1 } else { -3 });
//     let alns: Vec<_> = queries
//         .iter()
//         .map(|query| {
//             let aln = aligner.global(&query, &template);
//             let mut seq = vec![];
//             let mut ins = vec![b'-'; template.len()];
//             let (mut rpos, mut qpos) = (0, 0);
//             let mut prev = None;
//             for &op in aln.operations.iter() {
//                 use bio::alignment::AlignmentOperation::*;
//                 match op {
//                     Del => {
//                         seq.push(b'-');
//                         rpos += 1;
//                     }
//                     Ins => {
//                         if prev != Some(Ins) {
//                             ins[rpos] = query[qpos];
//                         }
//                         qpos += 1;
//                     }
//                     Subst | Match => {
//                         seq.push(query[qpos]);
//                         rpos += 1;
//                         qpos += 1;
//                     }
//                     _ => panic!(),
//                 }
//                 prev = Some(op);
//             }
//             let mut result = vec![];
//             for (x, y) in seq.into_iter().zip(ins) {
//                 result.push(x);
//                 result.push(y);
//             }
//             result
//         })
//         .collect();
//     let len = alns.iter().map(|s| s.len()).max().unwrap();
//     let summary: Vec<_> = {
//         let mut count = vec![vec![0; 5]; len];
//         for s in alns.iter() {
//             for (i, x) in s.iter().enumerate() {
//                 let idx = match *x {
//                     b'A' => 0,
//                     b'C' => 1,
//                     b'G' => 2,
//                     b'T' => 3,
//                     b'-' => 4,
//                     _ => 4,
//                 };
//                 count[i][idx] += 1;
//             }
//         }
//         count
//             .into_iter()
//             .map(|mut cs| {
//                 let (argmax, max) = cs.iter().enumerate().max_by_key(|x| x.1).unwrap();
//                 let max = *max;
//                 cs.remove(argmax);
//                 let (next, next_max) = cs.iter().enumerate().max_by_key(|x| x.1).unwrap().clone();
//                 let max_base = match argmax {
//                     0 => b'A',
//                     1 => b'C',
//                     2 => b'G',
//                     3 => b'T',
//                     4 => b'-',
//                     _ => b'-',
//                 };
//                 let next_max_base = match next {
//                     0 => b'A',
//                     1 => b'C',
//                     2 => b'G',
//                     3 => b'T',
//                     4 => b'-',
//                     _ => b'-',
//                 };
//                 eprintln!(
//                     "{}\t{}\t{}\t{}",
//                     max_base as char, max, next_max_base as char, next_max
//                 );
//                 (max_base, next_max_base)
//             })
//             .collect()
//     };
//     let serialized: Vec<_> = alns
//         .iter()
//         .map(|s| {
//             assert_eq!(s.len(), len);
//             // let component = s.iter().flat_map(|x| {
//             //     let mut res = vec![0.; 5];
//             //     let idx = match *x {
//             //         b'A' => 0,
//             //         b'C' => 1,
//             //         b'G' => 2,
//             //         b'T' => 3,
//             //         b'-' => 4,
//             //         _ => 4,
//             //     };
//             //     res[idx] = 1.;
//             //     res
//             // });
//             let component = s.iter().zip(summary.iter()).map(|(x, (m, n))| {
//                 if x == m {
//                     1.
//                 } else if x == n {
//                     -1.
//                 } else {
//                     0.
//                 }
//             });
//             DVector::from_iterator(len, component)
//             //DVector::from_iterator(5 * len, component)
//         })
//         .collect();
//     // Mean vector.
//     let mean = serialized
//         .iter()
//         //.fold(DVector::zeros(5 * len), |x, y| x + y)
//         .fold(DVector::zeros(len), |x, y| x + y)
//         / serialized.len() as f64;
//     let covariance = serialized
//         .iter()
//         .map(|x| (x - mean.clone()) * (x - mean.clone()).transpose())
//         .fold(DMatrix::zeros(len, len), |x, y| x + y)
//         //.fold(DMatrix::zeros(5 * len, 5 * len), |x, y| x + y)
//         / serialized.len() as f64;
//     let eigens = covariance.clone().symmetric_eigen();
//     let mut eigen_and_eigenvec: Vec<_> = eigens
//         .eigenvectors
//         .column_iter()
//         .zip(eigens.eigenvalues.iter())
//         .collect();
//     eigen_and_eigenvec.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
//     eigen_and_eigenvec.reverse();
//     for (_, val) in eigen_and_eigenvec.iter() {
//         eprintln!("{}", val);
//     }
//     let pca_vectors = &eigen_and_eigenvec[..3];
//     // let features: Vec<Vec<_>> = serialized
//     //     .iter()
//     //     .map(|x| pca_vectors.iter().map(|(v, _)| x.dot(&v)).collect())
//     //     .collect();
//     for (x, (n, _)) in serialized.iter().zip(reads) {
//         let mapped: Vec<_> = pca_vectors.iter().map(|(v, _)| x.dot(&v)).collect();
//         println!("{},{},{},{}", n, mapped[0], mapped[1], mapped[2]);
//     }
// }

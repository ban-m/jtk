// use bio::alignment::pairwise::Aligner;
// use nalgebra::*;
// use rand::Rng;
// use rand_xoshiro::Xoshiro256PlusPlus;
use definitions::*;
use rayon::prelude::*;
use serde_json;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let now = std::time::Instant::now();
    eprintln!("Begin");
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(std::io::BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    use std::io::BufRead;
    let alns: Vec<_> = std::fs::File::open(&args[2])
        .map(std::io::BufReader::new)?
        .lines()
        .filter_map(|line| line.ok())
        .filter_map(|line| bio_utils::lasttab::LastTAB::from_line(&line))
        .collect();
    let buckets = haplotyper::encode::distribute(&alns);
    eprintln!("Encoding:{}", (std::time::Instant::now() - now).as_secs());
    let read_name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    //Filter needed.
    let encoded_reads: Vec<_> = ds
        .raw_reads
        .par_iter()
        .filter_map(|read| {
            let alns = buckets.get(&read.name)?;
            haplotyper::encode::encode(&read, alns, &ds.selected_chunks)
        })
        .collect();
    eprintln!("Encoded:{}", (std::time::Instant::now() - now).as_secs());
    for encoded in encoded_reads.iter() {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for node in encoded.nodes.iter() {
            *counts.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        if counts.values().any(|&val| val > 1) {
            let units: Vec<_> = encoded
                .nodes
                .iter()
                .map(|n| format!("({},{})", n.unit, n.is_forward as u8))
                .collect();
            let id = encoded.id;
            let name = read_name[&id];
            println!("{}\t{}\t{}", id, name, units.join("-"));
        }
        for (&(unit_id, strand), _) in counts.iter().filter(|&(_, &val)| val > 1) {
            let chunk_len = 150;
            let unit = ds.selected_chunks.iter().find(|u| u.id == unit_id).unwrap();
            for node in encoded
                .nodes
                .iter()
                .filter(|n| n.unit == unit_id && strand == n.is_forward)
            {
                println!("{}\t{}", unit_id, node.is_forward);
                let (q, a, r) = node.recover(&unit);
                for ((q, a), r) in q
                    .chunks(chunk_len)
                    .zip(a.chunks(chunk_len))
                    .zip(r.chunks(chunk_len))
                {
                    println!("{}", String::from_utf8_lossy(q));
                    println!("{}", String::from_utf8_lossy(a));
                    println!("{}", String::from_utf8_lossy(r));
                }
            }
        }
    }
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

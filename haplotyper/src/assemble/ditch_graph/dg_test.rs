#[cfg(test)]
mod tests {
    #![allow(dead_code)]
    use super::super::*;
    use definitions::ReadType;
    use definitions::{Edge, Node};
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128StarStar;
    // Raw data, expected result.
    type DataPack = (Vec<EncodedRead>, Vec<Vec<u8>>);
    fn gen_node(unit: u64, cluster: u64, is_forward: bool, seq: Vec<u8>) -> Node {
        let mut node = Node::default();
        node.unit = unit;
        node.cluster = cluster;
        node.is_forward = is_forward;
        node.seq = seq.to_vec().into();
        node
    }
    fn gen_edge(from: u64, to: u64, offset: i64, label: Vec<u8>) -> Edge {
        Edge {
            from,
            to,
            offset,
            label: label.into(),
        }
    }
    //const BASES: &'static [char] = &['A', 'C', 'G', 'T'];
    const BASES: &'static [u8] = b"ACGT";
    use rand::seq::SliceRandom;
    fn revcmp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&x| match x {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect()
    }
    fn gen_mock1() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        let seq: Vec<_> = (0..30)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let mut read = EncodedRead::default();
        read.nodes.push(gen_node(0, 1, true, seq[..10].to_vec()));
        read.nodes.push(gen_node(1, 1, true, seq[11..26].to_vec()));
        read.nodes.push(gen_node(2, 1, false, revcmp(&seq[20..27])));
        read.nodes.push(gen_node(3, 1, true, seq[29..30].to_vec()));
        read.edges.push(gen_edge(0, 1, 1, seq[10..11].to_vec()));
        read.edges.push(gen_edge(1, 2, -6, vec![]));
        read.edges.push(gen_edge(2, 3, 2, seq[27..29].to_vec()));
        let reads = vec![read];
        (reads, vec![seq.to_vec()])
    }
    fn gen_mock_large() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        let seq_arms: Vec<Vec<_>> = (0..4)
            .map(|_| {
                (0..550)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let seq_center: Vec<_> = (0..100)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..5 {
            let seq: Vec<_> = seq_arms[0][unit * 110..unit * 110 + 100].to_vec();
            read0.nodes.push(gen_node(unit as u64, 0, true, seq));
            let seq: Vec<_> = seq_arms[1][unit * 110..unit * 110 + 100].to_vec();
            read1.nodes.push(gen_node(unit as u64, 1, true, seq));
            let lab: Vec<_> = seq_arms[0][unit * 110 + 100..(unit + 1) * 110].to_vec();
            read0
                .edges
                .push(gen_edge(unit as u64, unit as u64 + 1, 10, lab));
            let lab = seq_arms[1][unit * 110 + 100..(unit + 1) * 110].to_vec();
            read1
                .edges
                .push(gen_edge(unit as u64, unit as u64 + 1, 10, lab));
        }
        read0.nodes.push(gen_node(5, 0, true, seq_center.clone()));
        read1.nodes.push(gen_node(5, 0, true, seq_center.clone()));
        let (mut read2, mut read3) = (read.clone(), read.clone());
        read2.nodes.push(gen_node(5, 0, true, seq_center.clone()));
        read3.nodes.push(gen_node(5, 0, true, seq_center.clone()));
        for unit in 0..5 {
            let lab = seq_arms[2][unit * 110..unit * 110 + 10].to_vec();
            read2
                .edges
                .push(gen_edge(unit as u64 + 5, unit as u64 + 6, 10, lab));
            let lab = seq_arms[3][unit * 110..unit * 110 + 10].to_vec();
            read3
                .edges
                .push(gen_edge(unit as u64 + 5, unit as u64 + 6, 10, lab));
            let seq = seq_arms[2][10 + unit * 110..(unit + 1) * 110].to_vec();
            read2.nodes.push(gen_node(unit as u64 + 6, 0, true, seq));
            let seq = seq_arms[3][10 + unit * 110..(unit + 1) * 110].to_vec();
            read3.nodes.push(gen_node(unit as u64 + 6, 1, true, seq));
        }
        let answer: Vec<_> = seq_arms
            .into_iter()
            .chain(std::iter::once(seq_center))
            .map(|cs| cs.into_iter().collect::<Vec<_>>())
            .collect();
        (vec![read0, read1, read2, read3], answer)
    }
    fn gen_mock_large_2() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(11231312);
        let seq: Vec<Vec<_>> = (0..=8)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        let edg = Vec::new();
        for unit in 0..7 {
            let (subseq, dir): (_, _) = if [0, 1, 4, 5, 7].contains(&unit) {
                (seq[unit].clone(), true)
            } else {
                (revcmp(&seq[unit]), false)
            };
            let (cl0, cl1) = if [0, 1, 6, 7].contains(&unit) {
                (0, 0)
            } else {
                (0, 1)
            };
            let unit = unit as u64;
            read0.nodes.push(gen_node(unit, cl0, dir, subseq.clone()));
            read1.nodes.push(gen_node(unit, cl1, dir, subseq));
            read0.edges.push(gen_edge(unit, unit + 1, 0, edg.clone()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, edg.clone()));
        }
        read0.nodes.push(gen_node(7, 0, true, seq[7].clone()));
        read1.nodes.push(gen_node(7, 0, true, seq[7].clone()));
        let (read0, read1) = (read0, read1);
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (1..=7).rev() {
            let (subseq, dir): (Vec<_>, _) = if [0, 1, 4, 5, 7].contains(&unit) {
                (revcmp(&seq[unit]).clone(), false)
            } else {
                (seq[unit].clone(), true)
            };
            let (cl0, cl1) = if [0, 1, 6, 7].contains(&unit) {
                (0, 0)
            } else {
                (0, 1)
            };
            let unit = unit as u64;
            read2.nodes.push(gen_node(unit, cl0, dir, subseq.clone()));
            read3.nodes.push(gen_node(unit, cl1, dir, subseq));
            read2.edges.push(gen_edge(unit - 1, unit, 0, edg.clone()));
            read3.edges.push(gen_edge(unit - 1, unit, 0, edg.clone()));
        }
        let seq0: Vec<_> = revcmp(&seq[0]);
        read2.nodes.push(gen_node(0, 0, false, seq0.clone()));
        read3.nodes.push(gen_node(0, 0, false, seq0.clone()));
        let reads = vec![read0, read1, read2, read3];
        let middle = (2..=5).fold(Vec::new(), |mut x, y| {
            x.extend(seq[y].iter().copied());
            x
        });
        let answer: Vec<Vec<u8>> = vec![
            seq[0].iter().chain(seq[1].iter()).copied().collect(),
            middle.clone(),
            middle,
            seq[6].iter().chain(seq[7].iter()).copied().collect(),
        ];
        (reads, answer)
    }
    fn gen_circle() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(111);
        let seq: Vec<_> = (0..880)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let mut read = EncodedRead::default();
        for unit in 0..8 {
            let subseq: Vec<u8> = seq[110 * unit..110 * unit + 100].to_vec();
            read.nodes.push(gen_node(unit as u64, 0, true, subseq));
            let lab: Vec<u8> = seq[110 * unit + 100..110 * (unit + 1)].to_vec();
            read.edges
                .push(gen_edge(unit as u64, (unit as u64 + 1) % 7, 10, lab));
        }
        let subseq: Vec<_> = seq[..100].to_vec();
        read.nodes.push(gen_node(0, 0, true, subseq));
        (vec![read], vec![seq.into_iter().collect()])
    }
    fn gen_complex() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(11212);
        let units: Vec<Vec<_>> = (0..18)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let mut reads = vec![];
        let (mut read0, mut read1) = (EncodedRead::default(), EncodedRead::default());
        for unit in 8..10 {
            let subseq = units[unit].clone();
            let unit = unit as u64;
            read0.nodes.push(gen_node(unit, 0, true, subseq.clone()));
            read1.nodes.push(gen_node(unit, 1, true, subseq.clone()));
            read0.edges.push(gen_edge(unit, unit + 1, 0, Vec::new()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, Vec::new()));
        }
        let unit = 10u64;
        let subseq = units[unit as usize].clone();
        read0.nodes.push(gen_node(unit, 0, true, subseq.clone()));
        read1.nodes.push(gen_node(unit, 0, true, subseq.clone()));
        reads.push(read0);
        reads.push(read1);
        let mut read = EncodedRead::default();
        for &unit in &[10u64, 11u64] {
            let subseq = units[unit as usize].clone();
            read.nodes.push(gen_node(unit, 0, true, subseq));
            read.edges.push(gen_edge(unit, unit + 1, 0, Vec::new()));
        }
        read.nodes.push(gen_node(12, 0, true, units[12].clone()));
        reads.push(read);
        let mut read = EncodedRead::default();
        read.nodes.push(gen_node(10, 0, true, units[10].clone()));
        read.edges.push(gen_edge(unit as u64, 13, 0, Vec::new()));
        for &unit in &[13, 14, 15, 16] {
            read.nodes
                .push(gen_node(unit, 0, true, units[unit as usize].clone()));
            read.edges.push(gen_edge(unit, unit + 1, 0, Vec::new()));
        }
        read.nodes.push(gen_node(17, 0, true, units[17].clone()));
        read.edges.push(gen_edge(17, 11, 0, Vec::new()));
        read.nodes.push(gen_node(11, 0, true, units[11].clone()));
        reads.push(read);
        let answer = vec![
            units[8].iter().chain(units[9].iter()).copied().collect(),
            units[8].iter().chain(units[9].iter()).copied().collect(),
            units[10].clone(),
            units[11].iter().chain(units[12].iter()).copied().collect(),
            units[13..=17].iter().flatten().copied().collect(),
        ];
        (reads, answer)
    }
    fn gen_mock_complex() -> DataPack {
        let (mut circle_reads, mut circle) = gen_circle();
        let (complex_reads, complex) = gen_complex();
        circle_reads.extend(complex_reads);
        circle.extend(complex);
        (circle_reads, circle)
    }
    fn gen_remove_test_1() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let edge = Vec::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..3 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1) = if unit == 1 { (0, 1) } else { (0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq));
        }
        for unit in 0..2 {
            read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
        }
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (0..3).rev() {
            let seq: Vec<_> = units[unit].clone();
            let seq: Vec<_> = revcmp(&seq);
            let unit = unit as u64;
            let (cl2, cl3) = if unit == 1 { (0, 1) } else { (0, 0) };
            read2.nodes.push(gen_node(unit, cl2, false, seq.clone()));
            read3.nodes.push(gen_node(unit, cl3, false, seq));
        }
        for unit in (0..2).rev() {
            read2.edges.push(gen_edge(unit + 1, unit, 0, edge.clone()));
            read3.edges.push(gen_edge(unit + 1, unit, 0, edge.clone()));
        }
        let answer: Vec<_> = units.iter().flatten().copied().collect();
        (vec![read0, read1, read2, read3], vec![answer])
    }
    fn gen_by_units(units: &[Vec<u8>], read: &[(u64, u64, bool)]) -> EncodedRead {
        let edge = Vec::new();
        let mut eread = EncodedRead::default();
        for &(u, c, b) in read.iter() {
            let seq = units[u as usize].clone();
            if b {
                eread.nodes.push(gen_node(u, c, b, seq));
            } else {
                let seq: Vec<u8> = seq.clone();
                // let seq: String = seq.into_iter().collect();
                eread.nodes.push(gen_node(u, c, b, seq));
            }
        }
        for w in read.windows(2) {
            eread.edges.push(gen_edge(w[0].0, w[1].0, 0, edge.clone()));
        }
        eread
    }
    fn gen_remove_test_2() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<Vec<u8>> = (0..5)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let reads = vec![
            vec![(0, 0, true), (1, 0, true), (2, 0, true)],
            vec![(2, 0, false), (1, 0, false), (0, 0, false)],
            vec![(0, 0, true), (1, 1, true), (2, 0, true)],
            vec![(2, 0, false), (1, 1, false), (0, 0, false)],
            vec![(1, 1, true), (4, 0, true), (3, 0, true)],
            vec![(3, 0, false), (4, 0, false), (1, 1, false)],
        ];
        let reads: Vec<_> = reads
            .into_iter()
            .map(|r| gen_by_units(&units, &r))
            .collect();
        let answers = vec![
            units[0..2].iter().flatten().copied().collect(),
            units[2].clone(),
            units[4].iter().chain(units[3].iter()).copied().collect(),
        ];
        (reads, answers)
    }
    fn gen_remove_test_3() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<Vec<_>> = (0..3)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let edge = Vec::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1, mut read2) = (read.clone(), read.clone(), read.clone());
        for unit in 0..3 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1, cl2) = if unit == 1 { (0, 1, 2) } else { (0, 0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq.clone()));
            read2.nodes.push(gen_node(unit, cl2, true, seq));
            if unit < 2 {
                read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read2.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            }
        }
        let (mut read3, mut read4, mut read5) = (read.clone(), read.clone(), read.clone());
        for unit in (0..3).rev() {
            let seq: Vec<_> = units[unit].clone();
            let seq: Vec<_> = revcmp(&seq);
            let unit = unit as u64;
            let (cl3, cl4, cl5) = if unit == 1 { (0, 1, 2) } else { (0, 0, 0) };
            read3.nodes.push(gen_node(unit, cl3, false, seq.clone()));
            read4.nodes.push(gen_node(unit, cl4, false, seq.clone()));
            read5.nodes.push(gen_node(unit, cl5, false, seq.clone()));
            if 0 < unit {
                read3.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read4.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read5.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
            }
        }
        let reads = vec![read0, read1, read2, read3, read4, read5];
        let answers = units.into_iter().flatten().collect();
        (reads, vec![answers])
    }
    fn gen_remove_test_4() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<Vec<_>> = (0..5)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let reads = vec![
            vec![(0, 0, true), (1, 0, true), (2, 0, true)],
            vec![(0, 0, true), (1, 1, true), (2, 0, true)],
            vec![(0, 0, true), (3, 0, true), (4, 0, true)],
            vec![(4, 0, false), (3, 0, false), (0, 0, false)],
            vec![(2, 0, false), (1, 0, false), (0, 0, false)],
        ];
        let reads: Vec<_> = reads.iter().map(|r| gen_by_units(&units, r)).collect();
        let answers = vec![
            units[0].clone(),
            units[1].iter().chain(units[2].iter()).copied().collect(),
            units[3].iter().chain(units[4].iter()).copied().collect(),
        ];
        (reads, answers)
    }
    fn gen_remove_test_5() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<Vec<_>> = (0..5)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let edge = vec![];
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..5 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1) = if unit == 1 { (0, 1) } else { (0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq));
            if unit < 4 {
                read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            }
        }
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (0..5).rev() {
            let seq: Vec<_> = units[unit].clone();
            let seq = revcmp(&seq);
            let unit = unit as u64;
            let (cl2, cl3) = if unit == 1 { (0, 1) } else { (0, 0) };
            read2.nodes.push(gen_node(unit, cl2, false, seq.clone()));
            read3.nodes.push(gen_node(unit, cl3, false, seq));
            if 0 < unit {
                read2.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read3.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
            }
        }
        let answer = units.iter().flatten().copied().collect();
        (vec![read0, read1, read2, read3], vec![answer])
    }
    #[test]
    fn create() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let _graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
    }
    #[test]
    fn generate() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        let _ = graph.spell(&c);
    }
    #[test]
    fn validate() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        let (segment, _edge, _group, _) = graph.spell(&c);
        eprintln!("{:?}", segment);
        eprintln!("{:?}", answer);
        assert_eq!(segment.len(), 1);
        let seq = segment[0].sequence.as_ref().unwrap().as_bytes();
        let revcmp = revcmp_str(segment[0].sequence.as_ref().unwrap());
        let is_the_same = seq == answer[0] || answer[0] == revcmp.as_bytes();
        assert!(is_the_same);
    }
    #[test]
    fn validate_large() {
        let c = AssembleConfig::default();
        let (reads, _answer) = gen_mock_large();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 4);
        // for ans in answer {
        //     let forward = segments
        //         .iter()
        //         .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
        //     let rev_ans: Vec<_> = ans.chars().collect();
        //     let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
        //     let reverse = segments
        //         .iter()
        //         .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
        //     assert!(forward || reverse, "Couldn't find {}", ans);
        // }
    }
    #[test]
    fn validate_large_2() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_large_2();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 4);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(
                forward || reverse,
                "Couldn't find {}",
                std::str::from_utf8(&ans).unwrap()
            );
        }
    }
    #[test]
    fn validate_complex() {
        let c = AssembleConfig::default();
        let (reads, _answer) = gen_mock_complex();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 6);
        assert_eq!(segments.len(), 6);
    }
    #[test]
    fn validate_remove_1() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_1();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(
                forward || reverse,
                "Couldn't find {}",
                std::str::from_utf8(&ans).unwrap()
            );
        }
    }
    #[test]
    fn validate_remove_2() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_2();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 2);
        assert_eq!(segments.len(), 3);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(
                forward || reverse,
                "Couldn't find {}",
                std::str::from_utf8(&ans).unwrap()
            );
        }
    }
    #[test]
    fn validate_remove_3() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_3();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(
                forward || reverse,
                "Couldn't find {}",
                std::str::from_utf8(&ans).unwrap()
            );
        }
    }
    #[test]
    fn validate_remove_4() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_4();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 2);
        assert_eq!(segments.len(), 3);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(forward || reverse, "Couldn't fin");
        }
    }
    #[test]
    fn validate_remove_5() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_5();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, ReadType::ONT, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![].into()));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![].into()));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![].into());
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &ans);
            let rev_ans = revcmp(&ans);
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap().as_bytes() == &rev_ans);
            assert!(forward || reverse, "Couldn't find");
        }
        assert!(graph.sanity_check());
    }
}

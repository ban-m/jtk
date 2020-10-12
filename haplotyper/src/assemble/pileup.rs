#![allow(dead_code)]
use bio_utils::lasttab::LastTAB;
use std::collections::HashMap;

const BASE_TABLE: [usize; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

pub struct Pileups<'a> {
    template: &'a [u8],
    matches: Vec<[u32; 4]>,
    insertions: Vec<[u32; 4]>,
    deletions: Vec<u32>,
    length: usize,
}

use definitions::RawRead;
impl<'a> Pileups<'a> {
    pub fn convert_into_pileup(
        alignment: &[LastTAB],
        segment: &'a gfa::Segment,
        reads: &[&RawRead],
        _c: &super::AssembleConfig,
    ) -> Self {
        let template = segment.sequence.as_ref().map(|e| e.as_bytes()).unwrap();
        let length = template.len();
        let mut deletions = vec![0; length];
        let mut insertions = vec![[0; 4]; length + 1];
        let mut matches = vec![[0; 4]; length];
        for (idx, &b) in template.iter().enumerate() {
            matches[idx][BASE_TABLE[b as usize]] += 1;
        }
        let reads: HashMap<_, &RawRead> = reads.iter().map(|&r| (r.name.clone(), r)).collect();
        for aln in alignment {
            if let Some(read) = reads.get(aln.seq2_name()) {
                let seq = if aln.seq2_direction().is_forward() {
                    read.seq().to_vec()
                } else {
                    bio_utils::revcmp(read.seq())
                };
                let (mut rpos, mut qpos) = (aln.seq1_start(), aln.seq2_start());
                use bio_utils::lasttab::Op;
                for op in aln.alignment() {
                    match op {
                        Op::Match(l) => {
                            for (i, &b) in seq[qpos..qpos + l].iter().enumerate() {
                                matches[i + rpos][BASE_TABLE[b as usize]] += 1;
                            }
                            qpos += l;
                            rpos += l
                        }
                        Op::Seq1In(l) => {
                            for &b in seq[qpos..qpos + l].iter() {
                                insertions[rpos][BASE_TABLE[b as usize]] += 1;
                            }
                            qpos += l;
                        }
                        Op::Seq2In(l) => {
                            deletions[rpos] += l as u32;
                            rpos += l;
                        }
                    }
                }
            }
        }
        Self {
            length,
            template,
            matches,
            deletions,
            insertions,
        }
    }
    pub fn generate(&self) -> Vec<u8> {
        let mut pos = 0;
        let mut seq = vec![];
        while pos < self.length {
            let start = pos;
            let base = self.template[start];
            let end = start
                + self.template[pos..]
                    .iter()
                    .take_while(|&&b| b == base)
                    .count();
            let total_coverage = self.matches[start..end]
                .iter()
                .map(|x| x.iter().sum::<u32>())
                .sum::<u32>();
            let coverage = total_coverage / (end - start) as u32;
            // let base_count = self.matches[start..end]
            //     .iter()
            //     .map(|x| x[BASE_TABLE[base as usize]])
            //     .sum::<u32>()
            //     + self.insertions[start..=end]
            //         .iter()
            //         .map(|x| x[BASE_TABLE[base as usize]])
            //         .sum::<u32>();
            // for i in start..end {
            //     debug!(
            //         "DUMP\t{}\t{}\t{:?}\t{:?}\t{}",
            //         i, base as char, self.matches[i], self.insertions[i], self.deletions[i]
            //     );
            // }
            let ins = self.insertions[start..=end]
                .iter()
                .map(|x| x[BASE_TABLE[base as usize]])
                .sum::<u32>();
            let rep_num = if ins > coverage / 2 { 1 } else { 0 };
            //let rep_num = (base_count + coverage / 2) / coverage;
            // debug!(
            //     "DET\t{}\t{}\t{}\t{}\t{}\t{}",
            //     base_count, coverage, total_coverage, start, end, rep_num
            // );
            for _ in start..end {
                seq.push(base);
            }
            for _ in 0..rep_num {
                seq.push(base);
            }
            pos = end;
        }
        seq
    }
}

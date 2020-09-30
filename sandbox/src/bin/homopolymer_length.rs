use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reference = bio_utils::fasta::parse_into_vec(&args[1])?;
    let reads = bio_utils::fasta::parse_into_vec(&args[2])?;
    let alignment = bio_utils::sam::load_sam_file(&args[3])?;
    let reference = reference.iter().find(|r| r.id() == &args[4]).unwrap();
    let alignment: Vec<_> = alignment
        .into_iter()
        .filter(|aln| aln.r_name() == &args[4])
        .collect();
    let mut profile: Vec<HashMap<_, HashMap<_, usize>>> = vec![HashMap::new(); 4];
    let reads: HashMap<_, _> = reads.iter().map(|r| (r.id().to_string(), r)).collect();
    use rayon::prelude::*;
    let result: Vec<_> = alignment
        .iter()
        .filter(|aln| aln.is_primary() && aln.pos() > 0)
        .filter_map(|aln| {
            reads
                .get(aln.q_name())
                .map(|read| into_profile(read, aln, reference))
        })
        .collect();
    for aln in result {
        for (base, prf) in aln.into_iter().enumerate() {
            for (template, obss) in prf {
                for obs in obss {
                    *profile[base]
                        .entry(template)
                        .or_default()
                        .entry(obs)
                        .or_default() += 1;
                }
            }
        }
    }
    println!("Base\tAnswer\tObserve\tFraction");
    // A->C->G->T.
    for (idx, prof) in profile.iter().enumerate() {
        let base = match idx {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!(),
        };
        let mut prof: Vec<_> = prof.iter().collect();
        prof.sort_by_key(|x| x.0);
        for (temp, obss) in prof {
            let sum = obss.values().sum::<usize>();
            let mut obss: Vec<_> = obss.iter().collect();
            obss.sort_by_key(|x| x.0);
            for (obs, &count) in obss {
                let frac = count as f64 / sum as f64;
                println!("{}\t{}\t{}\t{}", base, temp, obs, frac);
            }
        }
    }
    Ok(())
}

use bio_utils::fasta::Record;
fn into_profile(
    read: &Record,
    aln: &bio_utils::sam::Sam,
    reference: &Record,
) -> Vec<HashMap<usize, Vec<usize>>> {
    let template = reference.seq();
    let length = template.len();
    let mut insertions = vec![[0; 4]; length + 1];
    let mut matches = vec![0; length];
    let seq = if aln.is_forward() {
        read.seq().to_vec()
    } else {
        bio_utils::revcmp(read.seq())
    };
    let (mut rpos, mut qpos) = (aln.pos() - 1, 0);
    // Record alignemnt
    use bio_utils::sam::Op;
    for op in aln.cigar() {
        match op {
            Op::Match(l) | Op::Mismatch(l) | Op::Align(l) => {
                for (i, &b) in seq[qpos..qpos + l].iter().enumerate() {
                    matches[i + rpos] = b;
                }
                qpos += l;
                rpos += l
            }
            Op::SoftClip(l) | Op::HardClip(l) => qpos += l,
            Op::Insertion(l) => {
                for &b in seq[qpos..qpos + l].iter() {
                    insertions[rpos][BASE_TABLE[b as usize]] += 1;
                }
                qpos += l;
            }
            Op::Deletion(l) => {
                rpos += l;
            }
            _ => {}
        }
    }
    let mut profile: Vec<HashMap<usize, Vec<usize>>> = vec![HashMap::new(); 4];
    let mut pos = aln.pos();
    while pos < rpos {
        let start = pos;
        let base = template[start];
        let end = start + template[pos..].iter().take_while(|&&b| b == base).count();
        let base_count = matches[start..end].iter().filter(|&&x| base == x).count()
            + insertions[start..=end]
                .iter()
                .map(|x| x[BASE_TABLE[base as usize]])
                .sum::<usize>();
        profile[BASE_TABLE[base as usize]]
            .entry(end - start)
            .or_default()
            .push(base_count);
        pos = end;
    }
    profile
}
const BASE_TABLE: [usize; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

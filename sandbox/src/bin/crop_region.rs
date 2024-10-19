use std::io::*;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut alignment_1 = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(bio_utils::sam::Sam::from_reader)?;
    alignment_1.records.retain(|aln| {
        let (s, e) = aln.refr_aligned_region();
        e - s > 10000
    });
    let mut alignment_2 = std::fs::File::open(&args[2])
        .map(BufReader::new)
        .map(bio_utils::sam::Sam::from_reader)?;
    alignment_2.records.retain(|aln| {
        let (s, e) = aln.refr_aligned_region();
        e - s > 10000
    });

    let cut_positions: Vec<_> = std::fs::File::open(&args[3])
        .map(BufReader::new)?
        .lines()
        .map_while(Result::ok)
        .map(|line| {
            let line: Vec<_> = line.split('\t').collect();
            let ctg_name = line[0].to_string();
            let start: usize = line[1].parse().unwrap();
            let end: usize = line[2].parse().unwrap();
            (ctg_name, start, end)
        })
        .collect();
    let refr_positions: Vec<_> = cut_positions
        .iter()
        .map(|(ctg_name, start, end)| {
            let (name_start, r_start) = alignment_1
                .records
                .iter()
                .find_map(|aln| get_corresponding_base(aln, ctg_name, *start))
                .unwrap();
            let (name_end, r_end) = alignment_1
                .records
                .iter()
                .find_map(|aln| get_corresponding_base(aln, ctg_name, *end))
                .unwrap();
            assert_eq!(name_start, name_end);
            (name_start, r_start, r_end)
        })
        .collect();
    assert!(refr_positions.iter().all(|c| c.0 == refr_positions[0].0));
    let ctg_name = refr_positions[0].0;
    let start = refr_positions.iter().map(|x| x.1).min().unwrap();
    let end = refr_positions.iter().map(|x| x.2).max().unwrap();
    let mut contig_names: Vec<_> = alignment_2.records.iter().map(|aln| aln.q_name()).collect();
    contig_names.sort();
    contig_names.dedup();
    for q_name in contig_names {
        // Reverse map.
        let (_, q_start) = alignment_2
            .records
            .iter()
            .filter(|aln| aln.q_name() == q_name)
            .find_map(|aln| find_region(aln, ctg_name, start))
            .unwrap();
        let (_, q_end) = alignment_2
            .records
            .iter()
            .filter(|aln| aln.q_name() == q_name)
            .find_map(|aln| find_region(aln, ctg_name, end))
            .unwrap();
        println!("{q_name}:{q_start}-{q_end}");
    }
    Ok(())
}

fn find_region<'a>(
    aln: &'a bio_utils::sam::Record,
    ref_name: &str,
    ref_target: usize,
) -> Option<(&'a str, usize)> {
    // Reference coordinate -> Query coordinate.
    let (ref_aln_start, ref_aln_end) = aln.refr_aligned_region();
    if ref_aln_start <= ref_target && ref_target < ref_aln_end && ref_name == aln.r_name() {
        let mut rpos = ref_aln_start;
        let mut qpos = 0;
        let mut qry_start = None;
        for op in aln.cigar() {
            match op {
                bio_utils::sam::Op::Align(l)
                | bio_utils::sam::Op::Match(l)
                | bio_utils::sam::Op::Mismatch(l) => {
                    if (rpos..rpos + l).contains(&ref_target) {
                        assert!(qry_start.is_none());
                        qry_start = Some(qpos + ref_target - rpos);
                    }
                    rpos += l;
                    qpos += l;
                }
                bio_utils::sam::Op::Deletion(l) => {
                    if (rpos..rpos + l).contains(&ref_target) {
                        assert!(qry_start.is_none());
                        qry_start = Some(qpos);
                    }
                    rpos += l;
                }
                bio_utils::sam::Op::SoftClip(l)
                | bio_utils::sam::Op::HardClip(l)
                | bio_utils::sam::Op::Insertion(l) => qpos += l,
                _ => {}
            }
        }
        let qry_name = aln.q_name();
        Some((qry_name, qry_start.unwrap()))
    } else {
        None
    }
}

fn get_corresponding_base<'a>(
    aln: &'a bio_utils::sam::Record,
    name: &str,
    query_target: usize,
) -> Option<(&'a str, usize)> {
    // Query Coordinate -> Reference Coordinate
    let (query_aln_start, query_aln_end) = aln.query_aligned_region();
    if aln.q_name() == name && query_aln_start <= query_target && query_target < query_aln_end {
        let mut rpos = aln.pos();
        let mut qpos = query_aln_start;
        let mut refr_target = None;
        for op in aln.cigar().iter() {
            match op {
                bio_utils::sam::Op::Match(l)
                | bio_utils::sam::Op::Mismatch(l)
                | bio_utils::sam::Op::Align(l) => {
                    if (qpos..qpos + l).contains(&query_target) {
                        refr_target = Some(rpos + query_target - qpos);
                    }
                    qpos += l;
                    rpos += l;
                }
                bio_utils::sam::Op::Deletion(l) => rpos += l,
                bio_utils::sam::Op::Insertion(l)
                | bio_utils::sam::Op::SoftClip(l)
                | bio_utils::sam::Op::HardClip(l) => {
                    if (qpos..qpos + l).contains(&query_target) {
                        refr_target = Some(rpos);
                    }
                    qpos += l;
                }
                _ => {}
            }
        }
        Some((aln.r_name(), refr_target.unwrap()))
    } else {
        None
    }
}

use definitions::*;
// use haplotyper::misc::cramers_v;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::squish_erroneous_clusters::*;
    let config = SquishConfig::default();
    ds.squish_erroneous_clusters(&config);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // let mut unit_pairs: HashMap<_, u32> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     for (i, n1) in read.nodes.iter().enumerate() {
    //         for n2 in read.nodes.iter().skip(i + 1) {
    //             let key = (n1.unit.min(n2.unit), n1.unit.max(n2.unit));
    //             *unit_pairs.entry(key).or_default() += 1;
    //         }
    //     }
    // }
    // let units: HashMap<_, _> = ds
    //     .selected_chunks
    //     .iter()
    //     .map(|n| (n.id, n.cluster_num))
    //     .collect();
    // unit_pairs.retain(|_, val| 5 < *val);
    // unit_pairs.retain(|(u1, u2), _| 1 < units[u1] && 1 < units[u2]);
    // use rayon::prelude::*;
    // let cramers_vs: Vec<_> = unit_pairs
    //     .par_iter()
    //     .filter_map(|(&(u1, u2), _)| check_correl(&ds, u1, u2).map(|x| (u1, u2, x)))
    //     .collect();
    // let mut unit_to_relvector: HashMap<_, Vec<_>> = HashMap::new();
    // for &(u1, u2, (rel, _)) in cramers_vs.iter() {
    //     unit_to_relvector.entry(u1).or_default().push((u2, rel));
    //     unit_to_relvector.entry(u2).or_default().push((u1, rel));
    // }
    // const THR: f64 = 0.5;
    // let stiff_units: HashSet<_> = unit_to_relvector
    //     .iter()
    //     .filter(|(_, rels)| rels.iter().any(|&(_, rel)| THR < rel))
    //     .map(|x| x.0)
    //     .collect();
    // let suspic_units: HashSet<_> = unit_to_relvector
    //     .iter()
    //     .filter(|(u, _)| !stiff_units.contains(u))
    //     .filter(|(_, rels)| {
    //         rels.iter()
    //             .any(|(to, rel)| stiff_units.contains(to) && *rel < THR)
    //     })
    //     .map(|x| x.0)
    //     .collect();
    // println!("ID\tType\tScore");
    // for unit in ds.selected_chunks.iter() {
    //     let id = unit.id;
    //     let assigned = match (stiff_units.contains(&id), suspic_units.contains(&id)) {
    //         (true, false) => "stiff",
    //         (false, true) => "suspic",
    //         (false, false) => "isolated",
    //         _ => panic!("{}", id),
    //     };
    //     println!("{id}\t{assigned}\t{}", unit.score);
    // }
    // println!("REL\tUnit1\tUnit2\tRel\tCount");
    // for (unit1, unit2, (rel, purged)) in cramers_vs.iter() {
    //     println!("REL\t{unit1}\t{unit2}\t{rel:.3}\t{purged}");
    // }
    // let scores: HashMap<_, _> = ds.selected_chunks.iter().map(|n| (n.id, n.score)).collect();
    // let mut ave_rel: HashMap<_, (f64, usize)> = HashMap::new();
    // for &(unit1, unit2, (rel, purged)) in cramers_vs.iter() {
    //     let slot = ave_rel.entry(unit1).or_default();
    //     slot.0 += rel * purged as f64;
    //     slot.1 += purged;
    //     let slot = ave_rel.entry(unit2).or_default();
    //     slot.0 += rel * purged as f64;
    //     slot.1 += purged;
    // }
    // println!("AVE\tUnit\tTotal\tAve\tScore");
    // for (unit, (sum, total)) in ave_rel {
    //     let score = scores[&unit];
    //     let ave = sum / total as f64;
    //     println!("AVE\t{unit}\t{total}\t{ave}\t{score}");
    // }
    Ok(())
}

// fn check_correl(ds: &DataSet, unit1: u64, unit2: u64) -> Option<(f64, usize)> {
//     let mut occs = vec![];
//     for read in ds.encoded_reads.iter() {
//         for node1 in read.nodes.iter().filter(|n| n.unit == unit1) {
//             for node2 in read.nodes.iter().filter(|n| n.unit == unit2) {
//                 occs.push((node1, node2));
//             }
//         }
//     }
//     corel(&occs)
// }

// fn corel(pairs: &[(&Node, &Node)]) -> Option<(f64, usize)> {
//     let occs: Vec<_> = pairs
//         .iter()
//         .map(|(n1, n2)| (n1.cluster as u32, n2.cluster as u32))
//         .collect();
//     let (first_slot_len, second_slot_len) = occs
//         .iter()
//         .fold((0, 0), |(fst, snd), &(x, y)| (fst.max(x), snd.max(y)));
//     if first_slot_len == 0 || second_slot_len == 0 {
//         None
//     } else {
//         Some((cramers_v(&occs), pairs.len()))
//     }
// }

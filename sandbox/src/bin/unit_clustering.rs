use definitions::*;
#[macro_use]
extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let unit: u64 = args[2].parse().unwrap();
    debug!("Started");
    let mut dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let c = haplotyper::ClusteringConfig::clr(&dataset, 2, 100, 1000, 10, false);
    let ref_unit = dataset
        .selected_chunks
        .iter()
        .find(|u| u.id == unit)
        .unwrap()
        .clone();
    use std::collections::HashMap;
    let id2name: HashMap<_, _> = dataset
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    let mut units: Vec<_> = dataset
        .encoded_reads
        .iter_mut()
        .filter_map(|read| {
            let id = read.id;
            read.nodes
                .iter_mut()
                .enumerate()
                .find(|(_, n)| n.unit == unit)
                .map(|(idx, n)| (id, idx, n))
        })
        .enumerate()
        // .inspect(|(idx, (id, _, _))| {
        //     let cl = if id2name[id].contains("hapA") { 1 } else { 0 };
        //     eprintln!("{}\t{}\t{}", idx, id2name[id], cl);
        // })
        .map(|(_, x)| x)
        .collect();
    haplotyper::unit_clustering(&mut units, &c, &ref_unit, 10);
    let mut result: Vec<_> = units
        .iter()
        .map(|(id, _, unit)| (&id2name[&id], unit.cluster))
        // .inspect(|(id, cl)| eprintln!("{}\t{}", id, cl))
        .collect();
    result.sort_by_key(|x| x.1);
    let mut counts: HashMap<_, usize> = HashMap::new();
    for (i, x) in result {
        if i.contains("hapA") {
            *counts.entry((x, 0)).or_default() += 1;
        } else {
            *counts.entry((x, 1)).or_default() += 1;
        }
    }
    for pred in 0..=1 {
        for ans in 0..=1 {
            print!("{}\t", counts.get(&(pred, ans)).unwrap_or(&0));
        }
    }
    println!();
    Ok(())
}

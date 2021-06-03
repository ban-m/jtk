use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let unit_id = 226;
    //let unit_id = 1;
    let ref_unit = ds.selected_chunks.iter().find(|u| u.id == unit_id).unwrap();
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
        .collect();
    let k = ref_unit.cluster_num;
    use std::collections::HashMap;
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.name.to_string()))
        .collect();
    use haplotyper::em_correction::*;
    let repeat_num = 5;
    let (new_clustering, lk, _cluster_num) = (0..repeat_num as u64)
        .map(|s| {
            let config = Config::new(repeat_num, unit_id * s, k, unit_id, 5);
            em_clustering(&reads, &config)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    eprintln!("LK:{:.3}", lk);
    let mut units = vec![String::new(); 9];
    units[4] = format!("{}", unit_id);
    for (read, (_, _, cluster)) in reads.iter().zip(new_clustering.iter()) {
        let mut context = vec![String::new(); 9];
        let index = read.nodes.iter().position(|n| n.unit == unit_id).unwrap();
        //context[4] = format!("{}", read.nodes[index].cluster);
        context[4] = format!("{}", cluster);
        let is_forward = read.nodes[index].is_forward;
        for (idx, node) in read.nodes.iter().skip(index + 1).enumerate().take(4) {
            let position = if is_forward { 5 + idx } else { 3 - idx };
            context[position] = format!("{}", node.cluster);
            units[position] = format!("{}", node.unit);
        }
        for (idx, node) in read.nodes.iter().take(index).rev().enumerate().take(4) {
            let position = if is_forward { 3 - idx } else { 5 + idx };
            context[position] = format!("{}", node.cluster);
            units[position] = format!("{}", node.unit);
        }
        let is_hapa = id2desc[&read.id].starts_with("hapA") as u8;
        println!("{}\t{}\t{}", is_hapa, read.id, context.join("\t"));
    }
    println!("U\t{}\t{}", 0, units.join("\t"));
    Ok(())
}

use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use std::collections::HashMap;
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        // .map(|read| (read.id, read.name.to_string()))
        .map(|read| (read.id, read.desc.to_string()))
        .collect();
    let unit_id: u64 = args[2].parse().expect("input unit id as 2nd argument.");
    let mut units = vec![String::new(); 9];
    units[4] = format!("{}", unit_id);
    for read in ds
        .encoded_reads
        .iter()
        .filter(|read| read.nodes.iter().any(|n| n.unit == unit_id))
    {
        let mut context = vec!["-".to_string(); 9];
        let index = read.nodes.iter().position(|n| n.unit == unit_id).unwrap();
        context[4] = format!("{}", read.nodes[index].cluster);
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
        // let is_hapa = id2desc[&read.id].starts_with("hapA") as u8;
        let is_hapa = id2desc[&read.id].contains("252v2") as u8;
        println!("{}\t{}\t{}", is_hapa, read.id, context.join("\t"));
    }
    println!("U\t{}\t{}", 0, units.join("\t"));
    Ok(())
}

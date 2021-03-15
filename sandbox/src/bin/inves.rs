use definitions::*;
use serde_json;
// use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let mut ds = haplotyper::encode::deletion_fill::correct_unit_deletion(ds);
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| definitions::Assignment::new(r.id, 0))
        .collect();
    use haplotyper::Assemble;
    let assemble_config = haplotyper::AssembleConfig::new(1, 100, false);
    eprintln!("Start assembling {} reads", ds.encoded_reads.len());
    eprintln!("Assembled reads.");
    let gfa = ds.assemble_as_gfa(&assemble_config);
    println!("{}", gfa);
    // let id2desc: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.desc.clone()))
    //     .collect();
    // for read in ds.encoded_reads.iter() {
    //     let line: Vec<_> = read.nodes.iter().map(|n| format!("{}", n.unit)).collect();
    //     println!("{}\t{}", id2desc[&read.id], line.join(":"));
    // }
    Ok(())
}

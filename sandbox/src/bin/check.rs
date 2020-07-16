use definitions::*;
#[macro_use]
extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("trace")).init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    use std::collections::HashMap;
    let id_to_asn: HashMap<_, _> = dataset
        .assignments
        .iter()
        .map(|asn| (asn.id, asn.cluster))
        .collect();
    let (unit, cluster, group) = (648, 2, 0);
    for read in dataset
        .encoded_reads
        .iter()
        .filter(|read| match id_to_asn.get(&read.id) {
            Some(&asn) if asn == group => read
                .nodes
                .iter()
                .any(|n| n.unit == unit && n.cluster == cluster),
            _ => false,
        })
    {
        let units: Vec<_> = read
            .nodes
            .iter()
            .map(|n| format!("{}:{}", n.unit, n.cluster))
            .collect();
        println!("{}\t{}", read.id, units.join(" "));
    }

    Ok(())
}

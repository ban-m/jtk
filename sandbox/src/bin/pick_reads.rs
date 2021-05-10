use definitions::*;
use std::collections::HashSet;
use std::io::BufReader;

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // let units: HashSet<u64> = std::fs::File::open(&args[2])
    //     .map(BufReader::new)?
    //     .lines()
    //     .filter_map(|x| x.ok())
    //     .filter_map(|x| {
    //         // let mut x = x.split('-');
    //         // let unit: u64 = x.next()?.parse().ok()?;
    //         // let cluster: u64 = x.next()?.parse().ok()?;
    //         // Some((unit, cluster))
    //         x.parse().ok()
    //     })
    //     .collect();
    let units: HashSet<u64> = vec![69, 148].into_iter().collect();
    let color = 40;
    use std::collections::HashMap;
    let id2name: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.desc.to_string()))
        .collect();
    for read in ds.encoded_reads.iter() {
        if read
            .nodes
            .iter()
            //.any(|n| units.contains(&(n.unit, n.cluster)))
            .any(|n| units.contains(&n.unit))
        {
            let line: Vec<_> = read
                .nodes
                .iter()
                .map(|n| {
                    // if units.contains(&(n.unit, n.cluster)) {
                    if units.contains(&n.unit) {
                        format!("\x1b[38;5;{}m{}\x1b[m", color, n.unit)
                    } else {
                        format!("{}", n.unit)
                    }
                })
                .collect();
            println!(
                "{}\t{}\t{}",
                id2name[&read.id],
                read.original_length,
                line.join(":")
            );
        }
    }
    Ok(())
}

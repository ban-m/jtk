use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let mut edge_count: HashMap<_, u32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for w in read.nodes.windows(2) {
            let from = w[0].unit;
            let to = w[1].unit;
            let tuple = if from < to {
                ((w[0].unit, w[0].cluster), (w[1].unit, w[1].cluster))
            } else {
                ((w[1].unit, w[1].cluster), (w[0].unit, w[0].cluster))
            };
            *edge_count.entry(tuple).or_default() += 1;
        }
    }
    for (((funit, fcl), (tunit, tcl)), count) in edge_count {
        println!("{}\t{}\t{}\t{}\t{}", funit, fcl, tunit, tcl, count);
    }
    Ok(())
}

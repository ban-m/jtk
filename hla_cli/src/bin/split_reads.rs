use bio_utils::fasta::{Record, Writer};
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
// Synopsis:
// ```bath
// cargo run --release --bin split_reads -- ${JSON} ${OUTPUT_DIR}
// ```
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let rdr = std::fs::File::open(&args[1]).map(BufReader::new)?;
    let ds: DataSet = serde_json::de::from_reader(rdr).unwrap();
    let output_dir = std::path::PathBuf::from(&args[2]);
    std::fs::create_dir_all(&output_dir)?;
    let id2cluster: HashMap<_, usize> = ds
        .assignments
        .iter()
        .map(|asn| (asn.id, asn.cluster))
        .collect();
    let mut writers: HashMap<usize, _> = id2cluster
        .values()
        .map(|&cluster| {
            let mut path = output_dir.clone();
            path.push(format!("{}.fasta", cluster));
            let file = std::fs::File::create(path).unwrap();
            let wtr = Writer::new(file);
            (cluster, wtr)
        })
        .collect();
    for read in ds.raw_reads.iter() {
        let desc = if read.desc.is_empty() {
            None
        } else {
            Some(read.desc.clone())
        };
        let record = Record::with_data(read.name.as_str(), &desc, read.seq.as_bytes());
        if let Some(cluster) = id2cluster.get(&read.id) {
            if let Some(wtr) = writers.get_mut(cluster) {
                wtr.write_record(&record)?;
            }
        }
    }
    Ok(())
}

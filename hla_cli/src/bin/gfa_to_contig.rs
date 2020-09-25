use bio_utils::fasta::{Record, Writer};
use gfa::*;
use std::io::BufReader;
// Synopsis:
// ```bath
// cargo run --release --bin gfa_to_contig -- ${GFA} ${OUTPUT_DIR}
// ```
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let rdr = std::fs::File::open(&args[1]).map(BufReader::new)?;
    let gfa = gfa::GFA::from_reader(rdr);
    let output_dir = std::path::PathBuf::from(&args[2]);
    std::fs::create_dir_all(&output_dir)?;
    let mut writers = std::collections::HashMap::new();
    for record in gfa.iter() {
        if let Content::Seg(segment) = &record.content {
            if let Some(seq) = &segment.sequence {
                let id = &segment.sid;
                let name = id.split('_').nth(1).unwrap();
                let name: usize = match name.parse() {
                    Err(why) => panic!("{:?}", why),
                    Ok(res) => res,
                };
                let record = Record::with_data(&id, &None, seq.as_bytes());
                writers
                    .entry(name)
                    .or_insert_with(|| {
                        let mut path = output_dir.clone();
                        path.push(format!("{}.contig.fasta", name));
                        println!("{:?}", path);
                        Writer::new(std::fs::File::create(&path).unwrap())
                    })
                    .write_record(&record)?;
            }
        }
    }
    Ok(())
}

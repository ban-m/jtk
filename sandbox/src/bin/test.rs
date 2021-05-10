use definitions::*;
// use std::collections::HashMap;
// use std::io::BufRead;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    for (i, read) in ds.encoded_reads.iter().enumerate() {
        let raw_read = ds.raw_reads.iter().find(|r| r.id == read.id).unwrap();
        let raw_seq: Vec<_> = raw_read
            .seq()
            .iter()
            .map(|x| x.to_ascii_uppercase())
            .collect();
        let mut seq = read.recover_raw_read();
        seq.iter_mut().for_each(|x| x.make_ascii_uppercase());
        for (xs, ys) in raw_seq.chunks(200).zip(seq.chunks(200)) {
            assert_eq!(String::from_utf8_lossy(xs), String::from_utf8_lossy(ys));
        }
    }
    Ok(())
}

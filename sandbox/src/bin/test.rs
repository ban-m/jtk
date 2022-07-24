use definitions::*;
use std::io::*;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use std::collections::HashMap;
    let raw_seq: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r)).collect();
    for read in ds.encoded_reads.iter() {
        let orig = read.recover_raw_read();
        assert!(orig.iter().all(u8::is_ascii_uppercase));
        let raw: Vec<_> = raw_seq[&read.id]
            .seq()
            .iter()
            .map(|x| x.to_ascii_uppercase())
            .collect();
        assert_eq!(orig, raw);
    }
    Ok(())
}

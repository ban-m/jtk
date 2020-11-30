use definitions::*;
// use std::collections::HashMap;
use std::io::BufReader;

fn main() -> std::io::Result<()> {
    let ds = get_input_file()?;
    // let args: Vec<_> = std::env::args().collect();
    // let unit: u64 = args[1].parse().unwrap();
    eprintln!("{}\t{}", ds.raw_reads.len(), ds.input_file);
    for read in ds.encoded_reads.iter()
    // .filter(|r| r.nodes.iter().any(|n| n.unit == unit))
    {
        let line: Vec<_> = read.nodes.iter().map(|n| format!("{}", n.unit)).collect();
        println!(
            "{}\t{}\t{}\t{}",
            line.join("-"),
            read.original_length,
            read.leading_gap.len(),
            read.trailing_gap.len()
        );
    }
    Ok(())
}

fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

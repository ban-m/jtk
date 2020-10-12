use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let start = std::time::Instant::now();
    let mut ds = get_input_file()?;
    let end = std::time::Instant::now();
    eprintln!("{:?}", end - start);
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|x| Assignment::new(x.id, 0))
        .collect();
    use haplotyper::assemble::*;
    let config = AssembleConfig::new(24, 100, false);
    let gfa = ds.assemble_as_gfa(&config);
    println!("{}", gfa);
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

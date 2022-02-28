fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let lengths: Vec<_> = bio_utils::fastq::parse_into_vec(&args[1])?
        .iter()
        .map(|r| r.len())
        .collect();
    let stdout = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(stdout.lock());
    use std::io::Write;
    for len in lengths {
        writeln!(&mut wtr, "{}\t{}", args[2], len)?;
    }
    Ok(())
}

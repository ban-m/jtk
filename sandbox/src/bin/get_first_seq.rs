fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let nth = bio::io::fasta::Reader::from_file(&args[1])?
        .records()
        .filter_map(|r| r.ok())
        .nth(0)
        .unwrap();
    println!(">{}\n{}", nth.id(), String::from_utf8_lossy(nth.seq()));
    Ok(())
}

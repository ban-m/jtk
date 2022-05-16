use std::collections::HashSet;
use std::io::BufRead;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ok_contigs: HashSet<_> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|line| bio_utils::paf::PAF::new(&line))
        .map(|paf| paf.qname)
        .collect();
    for line in std::fs::File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
    {
        match line.chars().next() {
            Some('S') => {
                let name = line.split('\t').nth(1).unwrap();
                if ok_contigs.contains(name) {
                    println!("{line}");
                }
            }
            Some('L') => {
                let mut fields = line.split('\t');
                let from = fields.nth(1).unwrap();
                let to = fields.nth(1).unwrap();
                if ok_contigs.contains(from) && ok_contigs.contains(to) {
                    println!("{line}");
                }
            }
            _ => {}
        }
    }
    Ok(())
}

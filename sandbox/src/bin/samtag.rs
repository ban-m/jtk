use std::io::BufRead;
use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let read_tag: HashMap<_, u64> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .map(|line| {
            let mut line = line.split('\t');
            let readid = line.next().unwrap().to_string();
            let hap: u64 = line.next().unwrap().parse().unwrap();
            (readid, hap)
        })
        .collect();
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    for line in reader.lines().filter_map(|line| line.ok()) {
        if line.starts_with('@') {
            println!("{line}");
        } else {
            print!("{line}");
            let id = line.split('\t').next().unwrap();
            match read_tag.get(id) {
                Some(hap) => println!("\tHP:i:{hap}"),
                None => println!(),
            }
        }
    }
    Ok(())
}

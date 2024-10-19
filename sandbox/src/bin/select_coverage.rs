use std::io::BufRead;

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::collections::HashMap;
    let target_regions: HashMap<_, _> = args[2]
        .split_whitespace()
        .map(|region| {
            if region.contains(':') {
                let (region, range) = region.rsplit_once(':').unwrap();
                let (start, end) = range.rsplit_once('-').unwrap();
                let start: usize = start.parse().unwrap();
                let end: usize = end.parse().unwrap();
                (region.to_string(), Some((start, end)))
            } else {
                (region.to_string(), None)
            }
        })
        .collect();
    for (name, range) in target_regions.iter() {
        eprintln!("{name},{range:?}");
    }
    for line in std::fs::File::open(&args[1])
        .map(std::io::BufReader::new)?
        .lines()
        .map_while(Result::ok)
    {
        if line.starts_with("contig") {
            println!("{line}");
            continue;
        }
        let mut fields = line.split_whitespace();
        let ctgname = fields.next().unwrap();
        let position: usize = fields.next().unwrap().parse().unwrap();
        let coverage: f64 = fields.next().unwrap().parse().unwrap();
        if let Some(range) = target_regions.get(ctgname) {
            match *range {
                None => println!("{line}"),
                Some((start, end)) if start <= position && position < end => {
                    println!("{ctgname}\t{}\t{coverage}", position - start);
                }
                _ => {}
            }
        }
    }
    Ok(())
}

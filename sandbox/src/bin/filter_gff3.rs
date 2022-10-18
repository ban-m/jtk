use std::io::*;
// cargo run --release --bin filter_gff3 -- <GFF3> <CHR> <START> <END> > <FILTERED_GFF3>
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("<GFF3> <CHR> <START> <END> <CHR_NAME>");
        return Ok(());
    }
    let chr = &args[2];
    let start: usize = args[3].parse().unwrap();
    let end: usize = args[4].parse().unwrap();
    let chr_output = &args[5];
    for line in std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
    {
        if line.starts_with('#') {
            println!("{line}");
        } else if line.starts_with(chr) {
            let mut line: Vec<_> = line.split('\t').collect();
            let beg: usize = line[3].parse().unwrap();
            let stop: usize = line[4].parse().unwrap();
            if line[0] == chr && start < beg && stop < end {
                // GFF3 is 1-index.
                let beg = format!("{}", beg - start + 1);
                let stop = format!("{}", stop - start + 1);
                line[0] = chr_output;
                line[3] = &beg;
                line[4] = &stop;
                println!("{}", line.join("\t"));
            }
        }
    }
    Ok(())
}

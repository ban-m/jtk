// Convert Longshot VCF to WhatsHap-compatible VCF.
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let file = std::fs::File::open(&args[1]).map(BufReader::new)?;
    let references = bio_utils::fasta::parse_into_vec(&args[2])?;
    let lines: Vec<_> = file.lines().filter_map(|x| x.ok()).collect();
    for (i, line) in lines.iter().enumerate() {
        if i == 2 {
            for refr in references.iter().filter(|refr| not_in(refr, &lines)) {
                let (id, length) = (refr.id(), refr.seq().len());
                println!("##contig=<ID={id},length={length}>");
            }
        }
        if line.starts_with('#') {
            println!("{}", line);
        } else {
            let mut fields: Vec<_> = line.split('\t').map(|x| x.to_string()).collect();
            if let Some(idx) = fields.iter().position(|field| field == "GT:GQ:PS:UG:UQ") {
                convert(&mut fields[idx + 1]);
            }
            println!("{}", fields.join("\t"));
        }
    }
    Ok(())
}
fn not_in(refr: &bio_utils::fasta::Record, lines: &[String]) -> bool {
    let contains = lines
        .iter()
        .filter(|line| line.starts_with("##contig"))
        .any(|line| {
            let line = line.trim_start_matches("##contig=<ID=");
            let line = line.split_once(',').unwrap().0;
            line.trim() == refr.id()
        });
    !contains
}

// Convert LongShot VCF-field into whatshap compatible one.
fn convert(field: &mut String) {
    let gq: u32 = field
        .split(':')
        .nth(1)
        .map(|x| x.parse::<f64>().unwrap().floor() as u32)
        .unwrap();
    let fields: Vec<_> = field
        .split(':')
        .enumerate()
        .map(|(i, x)| match i == 1 {
            true => format!("{}", gq),
            false => x.to_string(),
        })
        .collect();
    *field = fields.join(":");
}

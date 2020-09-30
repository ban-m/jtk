use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let sam: Vec<_> = bio_utils::sam::load_sam_file(&args[1]).unwrap();
    let reference: HashMap<_, _> = bio_utils::fasta::parse_into_vec(&args[2])
        .unwrap()
        .into_iter()
        .map(|r| (r.id().to_string(), r))
        .collect();
    for record in sam {
        let seq = record.seq().as_bytes();
        if seq.len() > 1 && record.pos() < 100 {
            let cigar = record.cigar();
            let refr = match reference.get(record.r_name()) {
                Some(res) => res.seq(),
                None => {
                    eprintln!("{}", record.r_name());
                    continue;
                }
            };
            let (q, o, r) = bio_utils::sam::recover_alignment(&cigar, seq, &refr, record.pos());
            for ((q, o), r) in q.chunks(150).zip(o.chunks(150)).zip(r.chunks(150)).take(1) {
                println!("{}", String::from_utf8_lossy(q));
                println!("{}", String::from_utf8_lossy(o));
                println!("{}\n", String::from_utf8_lossy(r));
            }
        }
    }
    Ok(())
}

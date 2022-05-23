use bio_utils::fasta::*;
fn main() -> std::io::Result<()> {
    // First argument would be the name of the binary...
    let mut records: Vec<_> = std::env::args()
        .skip(1)
        .flat_map(|file| parse_into_vec(file).unwrap())
        .collect();
    records.sort_by(|r1, r2| match r1.seq().len().cmp(&r2.seq().len()) {
        std::cmp::Ordering::Equal => r1.seq().cmp(r2.seq()),
        x => x,
    });
    let (mut prev_len, mut prev_seq): (usize, &[u8]) = (0, &[]);
    for record in records.iter() {
        if prev_len != record.seq().len() || prev_seq != record.seq() {
            println!("{record}");
            prev_len = record.seq().len();
            prev_seq = record.seq();
        }
    }
    Ok(())
}

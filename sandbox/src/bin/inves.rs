use kiley::sam::Sam;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let records = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(Sam::from_reader)?;
    let record = kiley::fasta::read_fasta(&Some(&args[2]))?;
    let (mut xpos, mut ypos) = (0, 0);
    let (xs, ys) = (&record[0].1, &record[1].1);
    let ofs = 5;
    for op in records.records[0].cigar().iter() {
        match op {
            kiley::sam::Op::Match(l) | kiley::sam::Op::Align(l) => {
                xpos += l;
                ypos += l;
            }
            kiley::sam::Op::Insertion(l) => {
                ypos += l;
                let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
                let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
                println!("{},{},{},I,{},{}", xpos, ypos, l, xseq, yseq);
            }
            kiley::sam::Op::Deletion(l) => {
                xpos += l;
                let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
                let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
                println!("{},{},{},D,{},{}", xpos, ypos, l, xseq, yseq);
            }
            kiley::sam::Op::Mismatch(l) => {
                xpos += l;
                ypos += l;
                let xseq = String::from_utf8_lossy(&xs[xpos - ofs..xpos + ofs]);
                let yseq = String::from_utf8_lossy(&ys[ypos - ofs..ypos + ofs]);
                println!("{},{},{},M,{},{}", xpos, ypos, l, xseq, yseq);
            }
            kiley::sam::Op::SoftClip(l) | kiley::sam::Op::HardClip(l) => {
                ypos += l;
            }
            _ => {}
        }
    }
    Ok(())
}

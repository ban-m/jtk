use std::io::*;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use kiley::sam::Sam;
    let records = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(Sam::from_reader)?;
    let sequences = kiley::fasta::read_fasta(&Some(&args[2]))?;
    fn parse_seqname(name: &str) -> Option<(&str, usize)> {
        let mut info = name.split(':');
        let seqname = info.next()?;
        let r_start: usize = info
            .next()?
            .split('-')
            .next()
            .and_then(|x| x.parse().ok())?;
        Some((seqname, r_start))
    }
    let ofs: usize = args.get(3).and_then(|x| x.parse().ok()).unwrap_or(3);
    println!("xpos,ypos,len,type,is_homopolymer,xseq,yseq");
    for record in records.records.iter() {
        let (rstart, rs, qstart, qs) = {
            let (rname, rpos) = parse_seqname(record.r_name()).unwrap();
            let (qname, qpos) = parse_seqname(record.q_name()).unwrap();
            let rs = &sequences.iter().find(|(id, _)| id == rname).unwrap().1;
            let qs = &sequences.iter().find(|(id, _)| id == qname).unwrap().1;
            let rpos = rpos + record.pos() - 1;
            (rpos, rs, qpos, qs)
        };
        let (mut rpos, mut qpos) = (rstart, qstart);
        for op in records.records[0].cigar().iter() {
            let rseq = String::from_utf8_lossy(&rs[rpos - ofs..rpos + ofs]);
            let qseq = String::from_utf8_lossy(&qs[qpos - ofs..qpos + ofs]);
            match op {
                kiley::sam::Op::Match(l) | kiley::sam::Op::Align(l) => {
                    rpos += l;
                    qpos += l;
                }
                kiley::sam::Op::Insertion(l) => {
                    let is_same = qs[qpos..qpos + l].iter().all(|&x| x == qs[qpos]);
                    let is_vague = qs[qpos - 1] == qs[qpos] || qs[qpos] == qs[qpos + l];
                    let is_hp = is_same && is_vague;
                    qpos += l;
                    print!("{},{},{},I,{},", rpos - rstart, qpos - qstart, l, is_hp,);
                    println!("{},{}", rseq, qseq);
                }
                kiley::sam::Op::Deletion(l) => {
                    let is_same = rs[rpos..rpos + l].iter().all(|&x| x == rs[rpos]);
                    let is_vague = rs[rpos - 1] == rs[rpos] || rs[rpos] == rs[rpos + l];
                    let is_hp = is_same && is_vague;
                    rpos += l;
                    print!("{},{},{},D,{},", rpos - rstart, qpos - qstart, l, is_hp);
                    println!("{},{}", rseq, qseq);
                }
                kiley::sam::Op::Mismatch(l) => {
                    qpos += l;
                    rpos += l;
                    print!("{},{},{},M,false,", rpos - rstart, qpos - qstart, l);
                    println!("{},{}", rseq, qseq);
                }
                kiley::sam::Op::SoftClip(l) | kiley::sam::Op::HardClip(l) => {
                    qpos += l;
                }
                _ => {}
            }
        }
    }
    Ok(())
}

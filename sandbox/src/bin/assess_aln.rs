use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let mut result: HashMap<_, (u64, u64, u64, u64)> = HashMap::new();
    let args: Vec<_> = std::env::args().collect();
    let sam = bio_utils::sam::load_sam_file(&args[1])?;
    use bio_utils::sam;
    let to_tuple = |&op: &sam::Op| match op {
        sam::Op::Match(l) | sam::Op::Align(l) => (l as u64, 0, 0, 0),
        sam::Op::Mismatch(l) => (0, l as u64, 0, 0),
        sam::Op::Insertion(l) => (0, 0, l as u64, 0),
        sam::Op::Deletion(l) => (0, 0, 0, l as u64),
        _ => {
            eprintln!("{:?}", op);
            (0, 0, 0, 0)
        }
    };
    let add = |(x, y, w, z), (a, b, c, d)| (x + a, y + b, w + c, z + d);
    for record in sam.iter() {
        let query = record.q_name();
        let refr = record.r_name();
        let (mat, mism, ins, del) = record.cigar().iter().map(to_tuple).fold((0, 0, 0, 0), add);
        if let Some(x) = result.get_mut(&(query, refr)) {
            x.0 += mat;
            x.1 += mism;
            x.2 += ins;
            x.3 += del;
        } else {
            result.insert((query, refr), (mat, mism, ins, del));
        }
    }
    for ((q_name, r_name), (mat, mism, ins, del)) in result {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            q_name, r_name, mat, mism, ins, del
        );
    }
    Ok(())
}

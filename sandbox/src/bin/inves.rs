use kiley::gphmm::*;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let seqs: Vec<_> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|line| line.ok())
        .map(|l| {
            l.split('\t')
                .nth(1)
                // .map(|s| kiley::padseq::PadSeq::new(s.as_bytes()))
                .unwrap()
                .to_string()
        })
        .collect();
    let template = seqs[0].as_bytes();
    let hmm = GPHMM::clr();
    profile_multi_deletion_banded_check(&hmm, template, seqs[1].as_bytes(), 100);
    // let (_, ops, _) = hmm.align(seqs[0].as_bytes(), seqs[1].as_bytes());
    // let (ta, oa, qa) = kiley::hmm::recover(seqs[0].as_bytes(), seqs[1].as_bytes(), &ops);
    // for ((ta, oa), qa) in ta.chunks(200).zip(oa.chunks(200)).zip(qa.chunks(200)) {
    //     println!("{}", String::from_utf8_lossy(ta));
    //     println!("{}", String::from_utf8_lossy(oa));
    //     println!("{}", String::from_utf8_lossy(qa));
    //     println!();
    // }
    // let prof = banded::ProfileBanded::new(&hmm, &seqs[0], &seqs[1], 100).unwrap();
    // let lk = prof.lk();
    // let mut deletions = prof.to_deletion_table(6);
    // deletions.iter_mut().for_each(|x| *x -= lk);
    // for (pos, lk) in deletions.iter().enumerate() {
    //     println!("{}\t{}\t{:.3}", pos / 5, pos % 5, lk);
    // }

    Ok(())
}

fn profile_multi_deletion_banded_check<T: HMMType>(
    model: &GPHMM<T>,
    xs: &[u8],
    ys: &[u8],
    radius: isize,
) {
    use banded::ProfileBanded;
    use kiley::padseq::PadSeq;
    let xs: Vec<_> = xs.iter().rev().copied().take(350).collect();
    let ys: Vec<_> = ys.iter().rev().copied().take(350).collect();
    let orig_xs: Vec<_> = xs.to_vec();
    let (xs, ys) = (PadSeq::new(xs), PadSeq::new(ys));
    let profile = ProfileBanded::new(model, &xs, &ys, radius).unwrap();
    let len = 6;
    let lk = profile.lk();
    let difftable = profile.to_deletion_table(len);
    println!("LK:{}", lk);
    println!(
        "{}",
        difftable.iter().fold(f64::NEG_INFINITY, |x, &y| x.max(y)) - lk
    );
    for (pos, diffs) in difftable.chunks(len - 1).enumerate() {
        let mut xs: Vec<_> = orig_xs.clone();
        xs.remove(pos);
        for (i, lkd) in diffs.iter().enumerate() {
            xs.remove(pos);
            let xs = PadSeq::new(xs.as_slice());
            let lk = model
                .likelihood_banded_inner(&xs, &ys, radius as usize)
                .unwrap();
            assert!((lk - lkd).abs() < 10f64, "{},{},{},{}", lk, lkd, pos, i);
        }
    }
}

#![allow(unused_imports)]
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let target = 1623;
    let unit = ds.selected_chunks.iter().find(|u| u.id == target).unwrap();
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let nodes: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| id2desc[&r.id].contains("000251v2"))
        .flat_map(|r| r.nodes.iter().filter(|n| n.unit == target))
        .collect();
    let (seqs, mut ops): (Vec<_>, Vec<_>) = nodes
        .iter()
        .map(|n| {
            let ops = haplotyper::local_clustering::ops_to_kiley_ops(&n.cigar);
            (n.seq(), ops)
        })
        .unzip();
    let hmm = haplotyper::local_clustering::get_tuned_model(&ds);
    let template = unit.seq();
    let take = seqs.len();
    let band = ds.read_type.band_width(template.len()) / 2;
    eprintln!("{band}");
    let cons = hmm.polish_until_converge_with_take(&template, &seqs, &mut ops, band, take);
    let ops = edlib_sys::global(&template, &cons);
    use kiley::Op;
    let etok = [Op::Match, Op::Ins, Op::Del, Op::Mismatch];
    let ops: Vec<_> = ops.iter().map(|&op| etok[op as usize]).collect();
    let (xr, ar, yr) = kiley::recover(&template, &cons, &ops);
    for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
        eprintln!("{}", String::from_utf8_lossy(xr));
        eprintln!("{}", String::from_utf8_lossy(ar));
        eprintln!("{}\n", String::from_utf8_lossy(yr));
    }
    let dist = ops.iter().filter(|&&x| x != kiley::Op::Match).count();
    eprintln!("CONS\t{dist}");
    // let (read_error_rate, unit_error_rate, sd) =
    //     haplotyper::encode::deletion_fill::estimate_error_rate(&ds, 0.15);
    // let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    // eprintln!("{}", sd);
    // let mut ok = 0;
    // println!("read\tunit\terror\texpected\tthr");
    // for read in ds.encoded_reads.iter() {
    //     for node in read.nodes.iter() {
    //         let (_, aln, _) = node.recover(units[&node.unit]);
    //         let diff = aln.iter().filter(|&&x| x != b'|').count();
    //         let error_rate = diff as f64 / aln.len() as f64;
    //         let expected = read_error_rate[read.id as usize] + unit_error_rate[node.unit as usize];
    //         let threshold = expected + 10f64 * sd;
    //         println!(
    //             "{}\t{}\t{}\t{}\t{}",
    //             read.id, node.unit, error_rate, expected, threshold
    //         );
    //         ok += (error_rate < threshold) as usize;
    //     }
    // }
    // let total: usize = ds.encoded_reads.iter().map(|r| r.nodes.len()).sum();
    // eprintln!("{}\t{}", ok, total);
    // for read in ds.encoded_reads.iter() {
    //     let error = read_error_rate[read.id as usize];
    //     let errors: f64 = read
    //         .nodes
    //         .iter()
    //         .map(|node| {
    //             let (_, aln, _) = node.recover(units[&node.unit]);
    //             let diff = aln.iter().filter(|&&x| x != b'|').count();
    //             diff as f64 / aln.len() as f64
    //         })
    //         .sum();
    //     let ave = errors / read.nodes.len() as f64;
    //     let len = read.nodes.len();
    //     let orig = read.original_length;
    //     println!("{}\t{}\t{}\t{}\t{}\tRead", read.id, error, ave, len, orig);
    // }
    // for chunk in ds.selected_chunks.iter() {
    //     let error = unit_error_rate[chunk.id as usize];
    //     println!("{}\t{}\tUnit", chunk.id, error);
    // }
    Ok(())
}

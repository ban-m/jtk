use definitions::*;
// use haplotyper::DetermineUnit;
use std::collections::HashMap;
use std::io::*;
fn error(node: &definitions::Node, ref_unit: &Unit) -> f64 {
    let (_, aln, _) = node.recover(ref_unit);
    let dist = aln.iter().filter(|&&x| x != b'|').count();
    dist as f64 / aln.len() as f64
    // let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
    // let del = query.iter().filter(|&&x| x == b' ').count() as f64;
    // let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
    // let aln_len = aln.len() as f64;
    // (mismat / aln_len, del / aln_len, ins / aln_len)
}
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    use rayon::prelude::*;
    let error_rates: Vec<_> = ds
        .encoded_reads
        .par_iter()
        .flat_map(|r| {
            r.nodes
                .par_iter()
                .map(|n| (r.id, n.unit, n.cluster, error(n, ref_units[&n.unit])))
                .collect::<Vec<_>>()
        })
        .collect();
    use haplotyper::determine_units::TAKE_THR;
    let sim_thr = haplotyper::determine_units::calc_sim_thr(&ds, TAKE_THR);
    // let sim_thr = haplotyper::determine_units::calc_sim_thr(&ds, 0.5);
    eprintln!("{sim_thr}");
    let (read_error_rate, unit_error_rate, median_of_var) =
        haplotyper::encode::deletion_fill::estimate_error_rate_dev(&ds, sim_thr);
    eprintln!("SD\t{median_of_var}");
    if let Ok(mut wtr) = std::fs::File::create("dump.log").map(BufWriter::new) {
        for read in ds.encoded_reads.iter() {
            let error = read_error_rate[read.id as usize];
            let len = read.original_length;
            writeln!(wtr, "READ\t{}\t{}\t{}", read.id, error, len)?;
        }
        for unit in ds.selected_chunks.iter() {
            for (cl, error) in unit_error_rate[unit.id as usize].iter().enumerate() {
                writeln!(wtr, "UNIT\t{}\t{cl}\t{error}", unit.id)?;
            }
        }
    }
    println!("ID\tUNIT\tCLUSTER\tERROR\tEXPECTED");
    for (id, unit, cluster, error) in error_rates {
        let expected =
            read_error_rate[id as usize] + unit_error_rate[unit as usize][cluster as usize];
        println!("{id}\t{unit}\t{cluster}\t{error}\t{expected}");
    }
    Ok(())
}

#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::ditch_graph::Focus;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use sandbox::IS_MOCK;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let min_span_read = ds.read_type.min_span_reads();
    let llr = ds.read_type.min_llr_value();
    let config = AssembleConfig::new(1, 1000, false, true, min_span_read, llr);
    check_foci(&ds, &config);
    Ok(())
}

pub fn check_foci(ds: &DataSet, c: &AssembleConfig) {
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, c);
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    let cov = ds.coverage.unwrap();
    graph.remove_zero_copy_elements(&lens, 0.2);
    graph.assign_copy_number(cov, &lens);
    graph.remove_zero_copy_elements(&lens, 0.5);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    let foci = graph.get_foci_dev(&reads, c);
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    // for focus in foci.iter().filter(|f| not_good(f, &counts)) {
    for focus in foci {
        let fcount = counts[&focus.from];
        let tcount = counts[&focus.to];
        println!("{focus}\t{:?}\t{:?}", fcount, tcount);
    }
}

// fn not_good(focus: &Focus, counts: &HashMap<(u64, u64), [u32; 2]>) -> bool {
//     let fcount = counts[&focus.from];
//     let tcount = counts[&focus.to];
//     (fcount[0] < fcount[1]) != (tcount[0] <= tcount[1])
// }

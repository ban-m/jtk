use definitions::DataSet;
use std::io::{BufRead, BufReader};
// const TARGET: &'static str = "hapD,N,684235-702082";
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let _lasttab: Vec<_> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|line| line.ok())
        .filter_map(|line| bio_utils::lasttab::LastTAB::from_line(&line))
        // .filter(|tab| tab.seq2_name() == TARGET)
        .collect();
    let _ds = get_input_file()?;
    // for read in ds.raw_reads.iter() {
    //     //let read = ds.raw_reads.iter().find(|r| r.name == TARGET).unwrap();
    //     let alignments: Vec<_> = lasttab
    //         .iter()
    //         .filter(|tab| tab.seq2_name() == read.name)
    //         .collect();
    // if let Some(result) = haplotyper::encode::encode(read, &alignments, &ds.selected_chunks) {
    // for node in result.nodes.iter() {
    //     for chunk in haplotyper::node_to_subchunks(&node, 100) {
    //         println!("{}\t{}", chunk.pos, chunk.seq.len());
    //     }
    // }
    // let len = result.nodes.iter().map(|n| n.seq.len()).sum::<usize>() as i64;
    // let edge_len = result.edges.iter().map(|n| n.offset).sum::<i64>();
    // let chunked_len = (len + edge_len) as usize;
    // println!("{}", result.original_length);
    // println!("{}", chunked_len + result.leading_gap + result.trailing_gap);
    // }
    // }
    Ok(())
}
fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

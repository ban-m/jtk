use definitions::*;
use std::io::BufReader;
fn main() {
    let ds = get_input_file().unwrap();
    // let max_unit = ds.selected_chunks.iter().map(|u| u.id).max().unwrap() as usize;
    // eprintln!("{}", max_unit);
    // let mut lengths: Vec<_> = (0..=max_unit).map(|x| vec![vec![]; x]).collect();
    // for read in ds.encoded_reads.iter() {
    //     for edge in read.edges.iter() {
    //         let (i, j) = (edge.from.max(edge.to), edge.from.min(edge.to));
    //         lengths[i as usize][j as usize].push(edge.offset);
    //     }
    // }
    // for (i, edges) in lengths.iter().enumerate() {
    //     for (j, ls) in edges.iter().enumerate() {
    //         if ls.len() >= 1 {
    //             let mean = ls.iter().sum::<i64>() / ls.len() as i64;
    //             eprintln!("{}\t{}\t{}\t{:.?}", i, j, mean, ls);
    //         }
    //     }
    // }
    for read in ds.encoded_reads {
        let line: Vec<_> = read
            .nodes
            .iter()
            .map(|n| format!("{}-{}", n.unit, n.cluster))
            .collect();
        println!("{}\t{}", read.id, line.join(","));
    }
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

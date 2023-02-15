#[allow(unused_macros)]
macro_rules! elapsed {
    ($a:expr) => {{
        let start = std::time::Instant::now();
        let return_value = $a;
        let end = std::time::Instant::now();
        (return_value, (end - start))
    }};
}

use definitions::*;

use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for read in ds.encoded_reads.iter_mut() {
        let seq = read.recover_raw_read();
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        fn start_end_iden(
            node: &definitions::Node,
            chunks: &HashMap<u64, &Chunk>,
        ) -> (usize, usize, f64) {
            let start = node.position_from_start;
            let chunk = chunks[&node.chunk];
            let (_, ar, _) = node.recover(chunk);
            let identity =
                1f64 - ar.iter().filter(|&&o| o != b'|').count() as f64 / ar.len() as f64;
            (start, start + node.seq().len(), identity)
        }
        loop {
            let should_be_removed = nodes.windows(2).enumerate().find_map(|(i, w)| {
                let (former_start, former_end, former_identity) = start_end_iden(&w[0], &chunks);
                let (latter_start, latter_end, latter_identity) = start_end_iden(&w[1], &chunks);
                if (former_start <= latter_start && latter_end < former_end)
                    || (latter_start <= former_start && former_end < latter_end)
                {
                    if former_identity < latter_identity {
                        Some(i)
                    } else {
                        Some(i + 1)
                    }
                } else {
                    None
                }
            });
            if let Some(idx) = should_be_removed {
                let removed = nodes.remove(idx);
                eprintln!("{}\t{idx}\t{removed}", read.id);
            } else {
                break;
            }
        }
        *read = haplotyper::encode::nodes_to_encoded_read(read.id, nodes, &seq).unwrap();
    }
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}

use definitions::*;
use std::collections::HashSet;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // let id2desc: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.desc.clone()))
    //     .collect();
    // let units: HashSet<(u64, u64)> = vec![(1317, 0)].into_iter().collect();
    let units: HashSet<u64> = vec![1317].into_iter().collect();
    for read in ds.encoded_reads.iter() {
        let mut dumps = vec![format!("{:<5}", 0); 7];
        for (idx, node) in read
            .nodes
            .iter()
            .enumerate()
            // .filter(|(_, n)| units.contains(&(n.unit, n.cluster)))
            .filter(|(_, n)| units.contains(&n.unit))
        {
            let (nodes, idx) = {
                let mut nodes: Vec<_> = read.nodes.iter().map(|n| n.unit).collect();
                match node.is_forward {
                    true => (nodes, idx),
                    false => {
                        nodes.reverse();
                        let idx = nodes.len() - idx - 1;
                        (nodes, idx)
                    }
                }
            };
            for i in (0..7).filter(|&i| 3 <= idx + i) {
                if let Some(unit) = nodes.get(idx + i - 3) {
                    dumps[i] = format!("{:<5}", unit);
                }
            }
            println!("{}\t{}", dumps.join("\t"), read.id);
        }
        // let pos: Vec<_> = id2desc[&read.id]
        //     .split(' ')
        //     .nth(0)
        //     .unwrap()
        //     .split(',')
        //     .collect();
        // let range: Vec<usize> = pos[2].split('-').map(|x| x.parse().unwrap()).collect();
        // let (start, end) = if pos[1] == "-strand" {
        //     (4_200_000 - range[1], 4_200_000 - range[0])
        // } else {
        //     (range[0], range[1])
        // };
        // if read.nodes.iter().any(|n| units.contains(&n.unit)) {
        // }
    }
    Ok(())
}

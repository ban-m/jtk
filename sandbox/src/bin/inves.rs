use definitions::*;

use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use std::collections::HashMap;
    let units: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    let num_cluster: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    println!("{units:?}");
    // collections::HashSet
    // let squished: HashSet<_> = ds
    //     .selected_chunks
    //     .iter()
    //     .filter_map(|c| (c.cluster_num == 1).then(|| c.id))
    //     .collect();
    for (i, &unit1) in units.iter().enumerate() {
        for &unit2 in units.iter().skip(i + 1) {
            let mut occs = vec![vec![0; num_cluster[&unit2]]; num_cluster[&unit1]];
            let (mut c1, mut c2) = (vec![], vec![]);
            for read in ds.encoded_reads.iter() {
                for node1 in read.nodes.iter().filter(|n| n.unit == unit1) {
                    for node2 in read.nodes.iter().filter(|n| n.unit == unit2) {
                        //if node1.is_biased(0.2) && node2.is_biased(0.2) {
                        occs[node1.cluster as usize][node2.cluster as usize] += 1;
                        c1.push(node1.cluster as usize);
                        c2.push(node2.cluster as usize);
                        // eprintln!("{:?}\t{:?}", node1.posterior, node2.posterior);
                        // let nodes: Vec<_> = read
                        //     .nodes
                        //     .iter()
                        //     .filter(|n| !squished.contains(&n.unit))
                        //     .map(|n| match n.is_biased(0.2) {
                        //         true => format!("{}-{}", n.unit, n.cluster),
                        //         false => format!("{}-?", n.unit),
                        //     })
                        //     .collect();
                        // eprintln!("{}", nodes.join("\t"));
                        //}
                    }
                }
                // let idx = read.nodes.iter().position(|n| n.unit == unit1);
                // let jdx = read.nodes.iter().position(|n| n.unit == unit2);
                // if let (Some(idx), Some(jdx)) = (idx, jdx) {
                //     eprintln!("{}", idx.max(jdx) - idx.min(jdx));
                // }
            }
            if occs.iter().flatten().sum::<u32>() == 0 {
                continue;
            }
            let ari = haplotyper::misc::adjusted_rand_index(&c1, &c2);
            println!("{unit1}\t{unit2}\t{ari:.3}");
            for row in occs.iter() {
                let row: Vec<_> = row.iter().map(|x| format!("{x}")).collect();
                println!("\t{}", row.join("\t"));
            }
        }
    }
    Ok(())
}

use std::{collections::HashMap, io::BufReader};

fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    ds.encoded_reads
        .iter_mut()
        .for_each(|r| r.nodes.iter_mut().for_each(|n| n.cluster = 0));
    let old_node_cp: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|x| (x.id, x.cluster_num))
        .collect();
    use haplotyper::multiplicity_estimation::MultiplicityEstimation;
    use haplotyper::multiplicity_estimation::MultiplicityEstimationConfig;
    let config = MultiplicityEstimationConfig::new(24, 4382094, Some("temp.gfa"));
    ds.estimate_multiplicity(&config);
    println!("NODE\tunit\tcluster\tcopy.number\tprev.estim");
    let mut node_count: HashMap<_, u32> = HashMap::new();
    let mut edge_count: HashMap<_, u32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *node_count.entry(node.unit).or_default() += 1
        }
        for edge in read.edges.iter() {
            *edge_count
                .entry((edge.from.min(edge.to), edge.from.max(edge.to)))
                .or_default() += 1;
        }
    }
    for (unit, cp) in ds.selected_chunks.iter().map(|c| (c.id, c.cluster_num)) {
        let old = old_node_cp[&unit];
        let count = node_count[&unit];
        println!("NODE\t{}\t{}\t{}\t{}", unit, count, cp, old);
    }
    // println!("EDGE\tfrom.unit\tfrom.cluster\tto.unit\tto.cluster\tcopy.number");
    // for (((fu, _), (tu, _)), cp) in edge_cp {
    //     let count = edge_count[&(fu.min(tu), fu.max(tu))];
    //     println!("EDGE\t{}\t{}\t{}\t{}", fu, tu, count, cp);
    // }
}

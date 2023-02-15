use definitions::*;
use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds_single: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let ds_both: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[2]).unwrap()))
            .unwrap();
    let mut cluster_num: HashMap<_, (_, _)> = HashMap::new();
    for c in ds_single.selected_chunks.iter() {
        cluster_num.entry(c.id).or_default().0 = c.cluster_num;
    }
    for c in ds_both.selected_chunks.iter() {
        cluster_num.entry(c.id).or_default().1 = c.cluster_num;
    }
    for (c, (single, both)) in cluster_num.iter() {
        println!("CLNUM\t{c}\t{single}\t{both}");
    }
    let mut clusterings: HashMap<_, (Vec<_>, Vec<_>)> = HashMap::new();
    for (r_single, r_both) in
        std::iter::zip(ds_single.encoded_reads.iter(), ds_both.encoded_reads.iter())
    {
        assert_eq!(r_single.nodes.len(), r_both.nodes.len());
        for (n_single, n_both) in std::iter::zip(r_single.nodes.iter(), r_both.nodes.iter()) {
            assert_eq!(n_single.chunk, n_both.chunk);
            let slot = clusterings.entry(n_single.chunk).or_default();
            slot.0.push(n_single.cluster as usize);
            slot.1.push(n_both.cluster as usize);
        }
    }
    for (uid, (c_single, c_both)) in clusterings {
        let rand_idx = haplotyper::misc::adjusted_rand_index(&c_single, &c_both);
        println!("RAND\t{uid}\t{rand_idx}");
    }
    Ok(())
}

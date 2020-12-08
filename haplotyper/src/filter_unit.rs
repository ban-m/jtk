use definitions::*;
// TODO: we should take parameter to trim tips, or short path.
pub trait FilterUnit {
    fn filter_unit(self, upper: usize, lower: usize) -> Self;
}

use std::collections::HashMap;
impl FilterUnit for DataSet {
    fn filter_unit(mut self, upper: usize, lower: usize) -> Self {
        let mut counts: HashMap<_, usize> = HashMap::new();
        for read in self.encoded_reads.iter() {
            for node in read.nodes.iter() {
                *counts.entry(node.unit).or_default() += 1;
            }
        }
        let original_len = self.selected_chunks.len();
        // debug!("REMOVED\tUnitID\tOcc");
        let filtered_unit: Vec<u64> = counts
            .iter()
            .filter_map(|(&key, &val)| {
                if lower < val && val < upper {
                    Some(key)
                } else {
                    // debug!("REMOVED\t{}\t{}", key, val);
                    None
                }
            })
            .collect();
        let filtered_unit: HashMap<u64, u64> = filtered_unit
            .into_iter()
            .enumerate()
            .map(|(idx, key)| (key, idx as u64))
            .collect();
        debug!("# of units:{}=>{}", original_len, filtered_unit.len());
        self.selected_chunks
            .retain(|chunk| filtered_unit.contains_key(&chunk.id));
        // Never panic.
        self.selected_chunks
            .iter_mut()
            .for_each(|chunk| chunk.id = *filtered_unit.get(&chunk.id).unwrap());
        let prev = self.encoded_reads.len();
        let node_prev = self
            .encoded_reads
            .iter()
            .map(|r| r.nodes.len())
            .sum::<usize>();
        self.encoded_reads.iter_mut().for_each(|read| {
            if log::log_enabled!(log::Level::Debug)
                && read
                    .nodes
                    .iter()
                    .all(|n| !filtered_unit.contains_key(&n.unit))
            {
                let read: Vec<_> = read
                    .nodes
                    .iter()
                    .map(|n| format!("{}({})", n.unit, counts.get(&n.unit).unwrap()))
                    .collect();
                debug!("{}", read.join("-"));
            }
            read.nodes
                .retain(|node| filtered_unit.contains_key(&node.unit));
            read.nodes
                .iter_mut()
                .for_each(|n| n.unit = *filtered_unit.get(&n.unit).unwrap());
        });
        self.encoded_reads.retain(|read| !read.nodes.is_empty());
        let now = self.encoded_reads.len();
        let node_now = self
            .encoded_reads
            .iter()
            .map(|r| r.nodes.len())
            .sum::<usize>();
        debug!(
            "Encoded Reads{}->{}(Raw reads{})",
            prev,
            now,
            self.raw_reads.len()
        );
        debug!("Number of Unit {}->{}", node_prev, node_now);
        self
    }
}

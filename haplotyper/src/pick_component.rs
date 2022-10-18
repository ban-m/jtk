use definitions::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct ComponentPickingConfig {
    /// How many component would we take?
    component_number: usize,
}

impl ComponentPickingConfig {
    pub fn new(component_number: usize) -> Self {
        Self { component_number }
    }
}

pub trait ComponentPicking {
    fn pick_top_n_component(&mut self, c: &ComponentPickingConfig);
}

impl ComponentPicking for DataSet {
    fn pick_top_n_component(&mut self, c: &ComponentPickingConfig) {
        use crate::assemble::AssembleConfig;
        // The last two parameter is not needed.
        let asm_config = AssembleConfig::new(100, false, false, 6, 3f64, false, None);
        let read_type = self.read_type;
        let reads = &self.encoded_reads;
        let units = &self.selected_chunks;
        let mut graph =
            crate::assemble::ditch_graph::DitchGraph::new(reads, units, read_type, &asm_config);
        debug!("CC\t{}", graph.cc());
        const LOWER_FRAC: f64 = 0.08;
        let cov = self.coverage.unwrap();
        let thr = (cov * LOWER_FRAC).round() as usize;
        graph.remove_lightweight_edges(thr, true);
        debug!("CC\t{}\tAfterRm", graph.cc());
        let mut components = graph.connected_components();
        for (i, cc) in components.iter().enumerate() {
            debug!("PICKING\t{}\t{}", i, cc.len());
        }
        components.sort_by_key(|cc| std::cmp::Reverse(cc.len()));
        let picked_units: HashSet<u64> = components
            .into_iter()
            .take(c.component_number)
            .flat_map(|cc| cc.iter().map(|&(unit, _)| unit).collect::<Vec<u64>>())
            .collect();

        debug!("PICKING\t{}\tConnectedComponents", c.component_number);
        debug!(
            "PICKING\t{}\t{}\tPicked",
            self.selected_chunks.len(),
            picked_units.len()
        );
        self.selected_chunks
            .retain(|unit| picked_units.contains(&unit.id));
        let mut unit_id_convert_table = HashMap::new();
        for unit in self.selected_chunks.iter_mut() {
            let new_id = unit_id_convert_table.len() as u64;
            unit_id_convert_table.insert(unit.id, new_id);
            unit.id = new_id;
        }
        let len = self.encoded_reads.len();
        self.encoded_reads.retain(|read| {
            read.nodes
                .iter()
                .all(|node| picked_units.contains(&node.unit))
        });
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                node.unit = unit_id_convert_table[&node.unit];
            }
            for edge in read.edges.iter_mut() {
                edge.from = unit_id_convert_table[&edge.from];
                edge.to = unit_id_convert_table[&edge.to];
            }
        }
        debug!(
            "PICKING\t{}\t{}\tEncodedRead",
            len,
            self.encoded_reads.len(),
        );
    }
}

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
    fn pick_top_n_component(self, c: &ComponentPickingConfig) -> Self;
}

impl ComponentPicking for DataSet {
    fn pick_top_n_component(mut self, c: &ComponentPickingConfig) -> Self {
        use crate::assemble::Assemble;
        use crate::assemble::AssembleConfig;
        let asm_config = AssembleConfig::new(1, 100, false, false, 6);
        let graph = self.assemble_draft_graph(&asm_config);
        let mut components: Vec<_> = graph.enumerate_connected_components();
        for (i, cc) in components.iter().enumerate() {
            let size = cc.iter().map(|node| node.segments.len()).sum::<usize>();
            debug!("PICKING\t{}\t{}", i, size);
        }
        components.sort_by_key(|cc| cc.iter().map(|node| node.segments.len()).sum::<usize>());
        components.reverse();
        let picked_units: HashSet<u64> = components
            .into_iter()
            .take(c.component_number)
            .flat_map(|cc| {
                cc.iter()
                    .flat_map(|node| node.segments.iter().map(|tile| tile.unit))
                    .collect::<Vec<u64>>()
            })
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
        let retained_reads: HashSet<_> = self.encoded_reads.iter().map(|r| r.id).collect();
        self.raw_reads
            .retain(|read| retained_reads.contains(&read.id));
        debug!(
            "PICKING\t{}\t{}\tEncodedRead",
            len,
            self.encoded_reads.len(),
        );
        // TODO: Maybe we need to fill the spase regions, right?
        self
    }
}

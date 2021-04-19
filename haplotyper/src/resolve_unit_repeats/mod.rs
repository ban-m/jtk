//! Resolve diploid level(**horizontal**) repeats by adjacency information.
//! For short, we collect all the adjacent units for each unit.
//! then, we construct a graph, where nodes are these adjacent units and edges are drawn
//! whenever these two node are connected by the picked unit (u -> unit -> v).
//! For each isolated node, we add edges from that isolated node to ALL OTHER node.
//! Then, for each connected component in that graph, we create new unit, with the same
//! sequence content, re-label all occurence of the node in the encoded reads.
use definitions::*;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct ResolveUnitRepeatsConfig {}

impl std::default::Default for ResolveUnitRepeatsConfig {
    fn default() -> Self {
        Self {}
    }
}

impl ResolveUnitRepeatsConfig {}

pub trait ResolveUnitRepeats {
    fn resolve_units(self, config: &ResolveUnitRepeatsConfig) -> Self;
}

impl ResolveUnitRepeats for DataSet {
    fn resolve_units(mut self, config: &ResolveUnitRepeatsConfig) -> Self {
        // Take unit id, return conversion information.
        let conversion_information: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|unit| get_conversion_information(unit, &self.encoded_reads, config))
            .collect();
        for read in self.encoded_reads.iter_mut() {
            let new_read_nodes = re_assign_units(&read, &conversion_information, config);
            read.nodes
                .iter_mut()
                .zip(new_read_nodes)
                .for_each(|(n, unit_id)| n.unit = unit_id);
        }
        let new_units: Vec<_> = conversion_information
            .values()
            .flat_map(|graph| graph.new_units())
            .collect();
        debug!("There would be {} new units.", new_units.len());
        self.selected_chunks.extend(new_units);
        self
    }
}

#[derive(Debug, Clone)]
struct AdjGraph {}

fn get_conversion_information(
    _unit: &Unit,
    _reads: &[EncodedRead],
    _c: &ResolveUnitRepeatsConfig,
) -> (u64, AdjGraph) {
    unimplemented!()
}

// Return the runs of the id of the units.
fn re_assign_units(
    _read: &EncodedRead,
    _infos: &HashMap<u64, AdjGraph>,
    _c: &ResolveUnitRepeatsConfig,
) -> Vec<u64> {
    unimplemented!()
}

impl AdjGraph {
    fn new_units(&self) -> Vec<Unit> {
        unimplemented!()
    }
}

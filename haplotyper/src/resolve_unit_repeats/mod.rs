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
    fn resolve_units(&mut self, config: &ResolveUnitRepeatsConfig);
}

impl ResolveUnitRepeats for DataSet {
    fn resolve_units(&mut self, config: &ResolveUnitRepeatsConfig) {
        // Take unit id, return conversion information.
        let mut current_max_unit = self.selected_chunks.iter().map(|x| x.id).max().unwrap();
        let conversion_information: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|unit| {
                let graph =
                    get_conversion_information(unit, &self.encoded_reads, current_max_unit, config);
                current_max_unit += graph.num_of_new_units();
                (unit.id, graph)
            })
            .collect();
        for read in self.encoded_reads.iter_mut() {
            re_assign_units(read, &conversion_information, config);
        }
        debug!(
            "ChunkNum\t{}\t{}",
            self.selected_chunks.len(),
            current_max_unit
        );
        for graph in conversion_information.values() {
            self.selected_chunks.extend(graph.new_units());
        }
    }
}

// Indicate how to convert each node based on adjacency information.
#[derive(Debug, Clone)]
struct AdjGraph {
    original_id: u64,
    // Map function.
    // It contains the original ID, also.
    // The second argument should be the direction of the adjacent node If you fix
    map: HashMap<(u64, bool), u64>,
}

fn get_conversion_information(
    _unit: &Unit,
    _reads: &[EncodedRead],
    _unit_id: u64,
    _c: &ResolveUnitRepeatsConfig,
) -> AdjGraph {
    unimplemented!()
}

// Return the runs of the id of the units.
fn re_assign_units(
    _read: &mut EncodedRead,
    _infos: &HashMap<u64, AdjGraph>,
    _c: &ResolveUnitRepeatsConfig,
) {
    unimplemented!()
}

impl AdjGraph {
    fn new_units(&self) -> Vec<Unit> {
        unimplemented!()
    }
    fn num_of_new_units(&self) -> u64 {
        unimplemented!()
    }
}

#![feature(is_sorted)]
pub mod assemble;
pub mod dense_encoding;
pub mod determine_units;
pub mod dirichlet_mixture;
pub mod encode;
pub mod entry;
pub mod extract;
pub mod find_union;
pub mod global_clustering;
pub mod hapcut;
pub mod local_clustering;
pub mod minimap2;
pub mod multiplicity_estimation;
pub mod pick_component;
pub mod polish_clustering;
pub mod polish_units;
pub mod re_clustering;
pub mod remove_erroneous_nodes;
pub mod repeat_masking;
pub mod stats;
pub mod unit_correlation;
pub mod view;
#[macro_use]
extern crate log;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

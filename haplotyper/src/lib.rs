#![feature(is_sorted)]
pub mod assemble;
pub mod consensus;
pub mod copy_number_estimation;
pub mod copy_number_estimation_mrf;
pub mod dense_encoding;
pub mod determine_units;
pub mod encode;
pub mod entry;
pub mod estimate_error_rate;
pub mod extract;
pub mod find_union;
pub mod likelihood_gains;
pub mod local_clustering;
pub mod minimap2;
pub mod misc;
pub mod model_tune;
pub mod multiplicity_estimation;
pub mod phmm_likelihood_correction;
pub mod pick_component;
pub mod polish_segments;
pub mod polish_units;
pub mod purge_diverged;
pub mod remove_erroneous_nodes;
pub mod repeat_masking;
pub mod seq;
pub mod stats;
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

pub const ALN_PARAMETER: (i32, i32, i32, i32) = (2, -6, -5, -1);

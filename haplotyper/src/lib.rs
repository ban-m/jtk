mod assemble;
mod determine_units;
pub mod em_correction;
pub mod encode;
mod entry;
mod extract;
mod find_union;
pub mod global_clustering;
pub mod local_clustering;
mod polish_clustering;
mod polish_units;
pub mod unit_correlation;
mod view;
#[macro_use]
extern crate log;
pub use assemble::*;
pub use encode::Encode;
pub use entry::Entry;
pub use extract::Extract;
pub use extract::ExtractTarget;
pub use view::View;
pub mod stats;
pub use determine_units::*;
pub use em_correction::ClusteringCorrection;
pub use global_clustering::*;
pub use local_clustering::*;
pub use polish_clustering::*;
pub use polish_units::*;
pub use stats::Stats;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

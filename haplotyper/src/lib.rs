mod assemble;
mod determine_units;
mod encode;
mod entry;
mod extract;
mod find_union;
mod global_clustering;
mod local_clustering;
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
pub use global_clustering::*;
pub use local_clustering::*;
pub use stats::Stats;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

mod determine_units;
mod encode;
mod extract;
pub use encode::Encode;
pub use extract::Extract;
pub use extract::ExtractTarget;
pub mod stats;
pub use stats::Stats;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

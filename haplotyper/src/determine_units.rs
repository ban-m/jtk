use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitConfig {}

pub trait DetermineUnit {
    fn select_chunks(self, config: UnitConfig) -> Self;
}

impl DetermineUnit for definitions::DataSet {
    fn select_chunks(self, config: UnitConfig) -> Self {
        self
    }
}

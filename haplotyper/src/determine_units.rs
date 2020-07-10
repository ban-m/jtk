use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitConfig {
    pub chunk_len: usize,
    pub chunk_num: usize,
    pub skip_len: usize,
    pub margin: usize,
}

pub trait DetermineUnit {
    fn select_chunks(self, config: &UnitConfig) -> Self;
}

use definitions::*;
impl DetermineUnit for definitions::DataSet {
    fn select_chunks(mut self, config: &UnitConfig) -> Self {
        let mut reads: Vec<&RawRead> = self.raw_reads.iter().collect();
        reads.sort_by_key(|r| r.seq().len());
        reads.reverse();
        debug!("Reads are sorted. {} reads.", reads.len());
        let selected_chunks: Vec<_> = reads
            .iter()
            .flat_map(|r| split_into(r, config))
            .take(config.chunk_num)
            .enumerate()
            .map(|(idx, c)| {
                let id = idx as u64;
                let seq = String::from_utf8_lossy(c).to_string();
                Unit { id, seq }
            })
            .collect();
        debug!(
            "Units collected. {}/{} units.",
            selected_chunks.len(),
            config.chunk_num
        );
        self.selected_chunks = selected_chunks;
        self
    }
}

fn split_into<'a>(r: &'a RawRead, c: &UnitConfig) -> Vec<&'a [u8]> {
    let seq = r.seq();
    if seq.len() < c.margin * 2 {
        vec![]
    } else {
        let end = seq.len() - c.margin;
        let stride = c.chunk_len + c.skip_len;
        (0..)
            .map(|i| (stride * i, stride * i + c.chunk_len))
            .map(|(x, y)| (x + c.margin, y + c.margin))
            .take_while(|&(_, y)| y < end)
            .map(|(s, t)| &seq[s..t])
            .collect()
    }
}

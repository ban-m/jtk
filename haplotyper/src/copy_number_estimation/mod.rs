pub mod copy_number_barrier;
pub mod copy_number_flow;
pub mod copy_number_gibbs;
pub mod copy_number_mrf;
pub mod copy_number_mst;

#[derive(Debug, Clone)]
pub struct CoverageCalibrator {
    // Sorted.
    lengths: Vec<usize>,
    // i->sum of the inverse value from i to the end.
    // cum_inverse_sum: Vec<f64>,
    // i->sum of the length of the reads longer than the i-th read.
    // Note that u64 is as big as 2^64 > total bp possibly.
    cum_sum: Vec<usize>,
    // Mean length.
    mean: f64,
}

impl CoverageCalibrator {
    pub fn new(lens: &[usize]) -> Self {
        let mut lengths = lens.to_vec();
        lengths.sort_unstable();
        let (mut cum_sum, sum): (Vec<_>, _) =
            lengths
                .iter()
                .rev()
                .fold((vec![], 0), |(mut sums, cum), &x| {
                    sums.push(cum + x);
                    (sums, cum + x)
                });
        cum_sum.reverse();
        let mean = sum as f64 / lengths.len() as f64;
        Self {
            lengths,
            cum_sum,
            mean,
        }
    }
    /// Calibrate the observed coverage in the gap-len bp region
    /// into the "actual" coverage, regarding that region as a single-bp region.
    /// Note that if the gap_len is longer than the longest reads in the dataset,
    /// it returns zero.
    pub fn calib(&self, observed: u64, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        let factor = self.cum_sum[idx] - gap_len * (self.lengths.len() - idx);
        observed as f64 * self.mean / (factor as f64 / self.lengths.len() as f64)
    }
    /// Calibrate the observed coverage in the gap-len bp region
    /// into the "actual" coverage, regarding that region as a single-bp region.
    /// Note that if the gap_len is longer than the longest reads in the dataset,
    /// it returns zero.
    pub fn calib_f64(&self, observed: f64, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        let factor = self.cum_sum[idx] - gap_len * (self.lengths.len() - idx);
        observed * self.mean / (factor as f64 / self.lengths.len() as f64)
    }
    /// Return the probability that a read span `gap_len` gap at a specific position.
    pub fn prob_spanning(&self, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        (self.cum_sum[idx] - gap_len * (self.lengths.len() - idx)) as f64
            / self.lengths.len() as f64
            / self.mean
    }
}

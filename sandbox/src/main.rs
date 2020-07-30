trait Count: std::marker::Sized {
    fn count(input: &[Self]) -> usize;
}

struct Test {
    inner: Vec<usize>,
}

impl Count for Test {
    fn count(input: &[Test]) -> usize {
        input.iter().map(|x| x.inner.len()).sum()
    }
}
fn main() {}

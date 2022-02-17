const fn revcmp() -> [u8; 256] {
    let mut table = [0; 256];
    table[b'A' as usize] = b'T';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'T' as usize] = b'A';
    table[b'a' as usize] = b't';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b't' as usize] = b'a';
    table
}
const REVCMP: [u8; 256] = revcmp();

pub struct DNAIter<'a> {
    index: usize,
    inner: &'a [u8],
    is_forward: bool,
}

impl<'a> DNAIter<'a> {
    pub fn new(seq: &'a [u8], is_forward: bool) -> Self {
        let mut it = Self {
            index: 0,
            inner: seq,
            is_forward,
        };
        if !is_forward {
            it.index = seq.len();
        }
        it
    }
}

impl<'a> std::iter::Iterator for DNAIter<'a> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.is_forward {
            self.index += 1;
            self.inner.get(self.index - 1).copied()
        } else if 0 < self.index {
            self.index -= 1;
            self.inner.get(self.index).map(|&x| REVCMP[x as usize])
        } else {
            None
        }
    }
}

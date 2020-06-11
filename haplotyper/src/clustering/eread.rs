pub struct ERead {
    pub id: u64,
    pub path: Vec<Elm>,
    pub cluster: usize,
}

pub struct Elm {
    pub unit: u64,
    pub cluster: usize,
}

impl Elm {
    pub fn new(unit: u64, cluster: usize) -> Self {
        Elm { unit, cluster }
    }
}

pub struct ChunkedUnit {
    pub cluster: usize,
    pub chunks: Vec<Chunk>,
}

pub struct Chunk {
    pub pos: usize,
    pub seq: Vec<u8>,
}

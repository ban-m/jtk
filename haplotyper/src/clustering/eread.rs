pub struct ERead {
    pub id: u64,
    pub path: Vec<Elm>,
    pub cluster: usize,
}

impl std::fmt::Debug for ERead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let path: Vec<_> = self
            .path
            .iter()
            .map(|p| format!("{}:{}", p.unit, p.cluster))
            .collect();
        write!(f, "{}\t{}\t{}", self.id, self.cluster, path.join("-"))
    }
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

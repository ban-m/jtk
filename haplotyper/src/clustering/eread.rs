pub struct ERead {
    pub id: u64,
    pub path: Vec<Elm>,
}

pub struct Elm {
    pub unit: u64,
    pub cluster: usize,
}

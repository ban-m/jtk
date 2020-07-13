use gfa::GFA;
pub trait Assemble {
    fn assemble_as_gfa(&self) -> GFA;
}

impl Assemble for definitions::DataSet {
    fn assemble_as_gfa(&self) -> GFA {
        unimplemented!()
    }
}

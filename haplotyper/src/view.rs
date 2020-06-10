pub trait View {
    fn view(&self, name: &str) -> Option<()>;
}

impl View for definitions::DataSet {
    fn view(&self, name: &str) -> Option<()> {
        let read = self.raw_reads.iter().find(|r| r.name == name)?;
        let encoded = self.encoded_reads.iter().find(|e| e.id == read.id)?;
        println!("Name:{}", name);
        println!("ID:{}", encoded.id,);
        println!("Length:{}", encoded.original_length);
        println!("{}-(Align)-{}", encoded.leading_gap, encoded.trailing_gap);
        println!("Encoded Rate:{:.3}", encoded.encoded_rate());
        for (idx, node) in encoded.nodes.iter().enumerate() {
            println!("{}-th unit, unit id is {}", idx, node.unit);
            println!("Start At {}", node.position_from_start);
            let d = if node.is_forward {
                "Forward"
            } else {
                "Reverse"
            };
            println!("Direction:{}", d);
            let unit = self.selected_chunks.iter().find(|u| u.id == node.unit)?;
            let (q, al, r) = node.recover(unit);
            for ((q, al), r) in q.chunks(50).zip(al.chunks(50)).zip(r.chunks(50)) {
                let q = String::from_utf8_lossy(q);
                let al = String::from_utf8_lossy(al);
                let r = String::from_utf8_lossy(r);
                println!("{}\n{}\n{}\n", q, al, r);
            }
        }
        for edge in encoded.edges.iter() {
            println!("{}-//{}//-{}", edge.from, edge.offset, edge.to);
        }
        println!("{}", read);
        Some(())
    }
}

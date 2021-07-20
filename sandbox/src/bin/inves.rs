use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    use haplotyper::encode::deletion_fill::ReadSkelton;
    let read_skeltons: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|read| read.nodes.len() > 1)
        .map(|read| ReadSkelton::from_rich_nodes(&read.nodes))
        .collect();
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.seq()))
        .collect();
    let selected_chunks = &ds.selected_chunks;
    let read_id = 8209;
    let r = ds
        .encoded_reads
        .iter_mut()
        .find(|r| r.id == read_id)
        .unwrap();
    let seq: Vec<_> = raw_seq[&r.id]
        .iter()
        .map(|x| x.to_ascii_uppercase())
        .collect();
    use haplotyper::encode::deletion_fill::correct_deletion_error;
    log::debug!("{}", r);
    correct_deletion_error(r, selected_chunks, &read_skeltons, &seq);
    // let target = 1374;
    // for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
    //     if edge.from == target || edge.to == target {
    //         let (x, y) = if edge.from == target {
    //             (edge.from, edge.to)
    //         } else {
    //             (edge.to, edge.from)
    //         };
    //         println!("{}\t{}\t{}", x, y, edge.offset);
    //     }
    // }
    // let unit = ds.selected_chunks.iter().find(|u| u.id == 375).unwrap();
    // for read in ds.encoded_reads.iter() {
    //     for (idx, node) in read.nodes.iter().enumerate().filter(|(_, n)| n.unit == 375) {
    //         let (_, ax, _) = node.recover(&unit);
    //         let dist = ax.iter().filter(|&&x| x != b'|').count();
    //         let to = match node.is_forward {
    //             true => read.nodes.get(idx + 1).map(|x| x.unit).unwrap_or(0),
    //             false => read.nodes.get(idx - 1).map(|x| x.unit).unwrap_or(0),
    //         };
    //         println!("DIST\t{}\t{}\t{}", unit.id, to, dist);
    // for ((qx, ax), rx) in qx.chunks(200).zip(ax.chunks(200)).zip(rx.chunks(200)) {
    //     let qx = String::from_utf8_lossy(qx);
    //     let ax = String::from_utf8_lossy(ax);
    //     let rx = String::from_utf8_lossy(rx);
    //     println!("{}\n{}\n{}\n", qx, ax, rx);
    // }
    //     }
    // }
    Ok(())
}

use definitions::*;
use sandbox::IS_MOCK;
use std::collections::HashMap;
use std::io::BufReader;
// Synopsis: [JSON] [Unit to dump]+
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    for query in &args[2..] {
        let mut query = query.split(':');
        let unit: u64 = query.next().unwrap().parse().unwrap();
        let range = query.next().map(|range| {
            let mut range = range.split('-');
            let start: usize = range.next().unwrap().parse().unwrap();
            let end: usize = range.next().unwrap().parse().unwrap();
            (start, end)
        });
        let ref_chunk = ds.selected_chunks.iter().find(|c| c.id == unit).unwrap();
        let mut nodes: Vec<_> = ds
            .encoded_reads
            .iter()
            .flat_map(|read| {
                let is_hap1 = match IS_MOCK {
                    true => id2desc[&read.id].contains("hapA") as usize,
                    false => id2desc[&read.id].contains("000251v2") as usize,
                };
                read.nodes
                    .iter()
                    .enumerate()
                    .filter(|n| unit == n.1.unit)
                    .map(|(i, n)| (read.id, is_hap1, i, n, read.nodes.len()))
                    .collect::<Vec<_>>()
            })
            .collect();
        nodes.sort_unstable_by_key(|x| (x.3.cluster, x.1));
        for (readid, is_hap1, i, node, len) in nodes {
            let (query, aln, refr) = node.recover(ref_chunk);
            let dist = aln.iter().filter(|&&x| x != b'|').count();
            let identity = 1f64 - dist as f64 / aln.len() as f64;
            let post: Vec<_> = node.posterior.iter().map(|p| format!("{:.2}", p)).collect();
            let (unit, cluster, post) = (node.unit, node.cluster, post.join("\t"));
            let name = &id2desc[&readid];
            println!(
                "{readid}\t{is_hap1}\t{i}\t{len}\t{unit}\t{cluster}\t{:.2}\t{post}\t{name}",
                identity,
            );
            println!("ALN\t{}", dist);
            let (start, end) = match range {
                Some((start, end)) => {
                    let (st_pos, _) = refr
                        .iter()
                        .enumerate()
                        .filter(|&(_, &b)| b != b' ')
                        .take(start)
                        .last()
                        .unwrap();
                    let (end_pos, _) = refr
                        .iter()
                        .enumerate()
                        .filter(|&(_, &b)| b != b' ')
                        .take(end)
                        .last()
                        .unwrap();
                    (st_pos, end_pos)
                }
                None => (0, refr.len()),
            };
            let alns = query.into_iter().zip(aln).zip(refr);
            let alns: Vec<_> = alns.skip(start).take(end - start).collect();
            for chunk in alns.chunks(200) {
                let query: Vec<_> = chunk.iter().map(|x| x.0 .0).collect();
                let aln: Vec<_> = chunk.iter().map(|x| x.0 .1).collect();
                let refr: Vec<_> = chunk.iter().map(|x| x.1).collect();
                println!("ALN\t{}", String::from_utf8_lossy(&query));
                println!("ALN\t{}", String::from_utf8_lossy(&aln));
                println!("ALN\t{}", String::from_utf8_lossy(&refr));
            }
        }
    }
    Ok(())
}

use definitions::*;
use sandbox::IS_MOCK;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use std::collections::HashMap;
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        for node in read.nodes.iter() {
            counts.entry((node.chunk, node.cluster)).or_default()[ans] += 1;
        }
    }
    let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
    println!("UNIT\tchunk\tcluster\thap1\thap2\tpurity\tscore");
    for ((chunk, cluster), counts) in counts.iter() {
        let score = score[chunk];
        let total = counts[0] + counts[1];
        let pur = counts[0].max(counts[1]) as f64 / total as f64;
        println!(
            "UNIT\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            chunk, cluster, counts[0], counts[1], pur, score,
        );
    }
    if args.len() < 3 {
        eprintln!("Contig information was not supplied. Finish after dumping chunk information.");
        return Ok(());
    }
    let tigs = std::fs::File::open(&args[2]).map(BufReader::new)?;
    use std::io::BufRead;
    let tigs: Vec<_> = tigs
        .lines()
        .map_while(Result::ok)
        .map(|line| {
            let mut line = line.split('\t');
            let id = line.next().unwrap().to_string();
            let copy_num: usize = line.next().unwrap().parse().unwrap();
            let clusters: Vec<_> = line
                .filter_map(|x| {
                    let mut elm = x.split('-');
                    let chunk: u64 = elm.next()?.parse().unwrap();
                    let cluster: u64 = elm.next()?.parse().unwrap();
                    Some((chunk, cluster))
                })
                .collect();
            (id, copy_num, clusters)
        })
        .collect();
    let mut max_cluster_id: HashMap<u64, u64> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        let _ = max_cluster_id
            .entry(node.chunk)
            .and_modify(|x| *x = (*x).max(node.cluster))
            .or_insert(node.cluster);
    }
    let counts: HashMap<String, (usize, [u32; 2])> = tigs
        .iter()
        .map(|&(ref id, copy_num, ref clusters)| {
            let cs = clusters
                .iter()
                .filter(|&(u, _)| copy_num != 1 || 0 < max_cluster_id[u])
                .filter_map(|key| counts.get(key))
                .fold([0, 0], |mut acc, xs| {
                    acc[0] += xs[0];
                    acc[1] += xs[1];
                    acc
                });
            (id.clone(), (copy_num, cs))
        })
        .collect();
    println!("CONTIG\tID\tCopyNum\tLen\tHap1\tHap2\tPurity");
    for (id, (copy_num, c)) in counts {
        let length = tigs.iter().find(|x| x.0 == id).map(|x| x.2.len()).unwrap();
        let total = c[0] + c[1];
        let pur = c[0].max(c[1]) as f64 / total as f64;
        println!(
            "CONTIG\t{}\t{}\t{}\t{}\t{}\t{:.4}",
            id, copy_num, length, c[0], c[1], pur,
        );
    }
    Ok(())
}

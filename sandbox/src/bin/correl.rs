use definitions::*;
use rayon::prelude::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let result: Vec<_> = ds
        .selected_chunks
        .par_iter()
        .flat_map(|unit| {
            ds.selected_chunks
                .iter()
                .filter_map(|jnit| {
                    let mut pairs = vec![];
                    let k1 = unit.cluster_num;
                    let k2 = jnit.cluster_num;
                    let mut count = 0;
                    for read in ds.encoded_reads.iter() {
                        let (in_i, in_j): (Vec<_>, Vec<_>) = read
                            .nodes
                            .iter()
                            .filter(|n| n.unit == unit.id || n.unit == jnit.id)
                            .partition(|n| n.unit == unit.id);
                        count += (!in_i.is_empty() && !in_j.is_empty()) as usize;
                        for inode in in_i {
                            for jnode in in_j.iter() {
                                pairs.push((inode.cluster, jnode.cluster));
                            }
                        }
                    }
                    (!pairs.is_empty() && 10 < count).then(|| {
                        let p_value = haplotyper::unit_correlation::calc_p_value(
                            &pairs, 4324, k1 as u64, k2 as u64,
                        );
                        (unit.id, jnit.id, k1, k2, p_value, count)
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect();

    println!("ID1\tID2\tCluste1\tCluster2\tCount\tPValue");
    for (id, jd, k1, k2, p_value, count) in result {
        println!("{}\t{}\t{}\t{}\t{}\t{}", id, jd, k1, k2, count, p_value);
    }
    Ok(())
}

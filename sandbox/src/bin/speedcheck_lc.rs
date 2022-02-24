const IS_MOCK: bool = true;
use definitions::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    println!("UNIT\tunit\tcluster\thap1\thap2\tpurity\tscore\ttype");
    println!("ELAPSE\tunit\ttime\tlen\tcoverage\tclusternum\ttype");
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .for_each(|n| n.cluster = 0);
    local_clustering(&mut ds, true);
    // Local clustering before
    dump(&ds, true);
    // Local clustering after
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .for_each(|n| n.cluster = 0);
    local_clustering(&mut ds, false);
    dump(&ds, false);
    Ok(())
}

fn dump(ds: &DataSet, key: bool) {
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
    for ((unit, cluster), counts) in counts.iter() {
        let score = score[unit];
        let total = counts[0] + counts[1];
        let pur = counts[0].max(counts[1]) as f64 / total as f64;
        println!(
            "UNIT\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
            unit, cluster, counts[0], counts[1], pur, score, key
        );
    }
}

fn local_clustering(ds: &mut DataSet, is_before: bool) {
    let mut pileups: HashMap<u64, Vec<&mut Node>> =
        ds.selected_chunks.iter().map(|c| (c.id, vec![])).collect();
    let chunks: HashMap<u64, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some(bucket) = pileups.get_mut(&node.unit) {
            bucket.push(node);
        }
    }
    pileups.iter_mut().for_each(|(unit_id, nodes)| {
        nodes.sort_by_cached_key(|node| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            let (_, aln, _) = node.recover(ref_unit);
            aln.iter().filter(|&&x| x != b'|').count()
        });
    });
    let coverage = ds.coverage.unwrap();
    let read_type = ds.read_type;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .filter(|(_, units)| !units.is_empty())
        .map(|(&unit_id, units)| {
            let ref_unit = chunks.get(&unit_id).unwrap();
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 25);
            let (seqs, mut ops): (Vec<_>, Vec<_>) = units
                .iter()
                .map(|node| {
                    let mut ops = vec![];
                    use definitions::Op;
                    for op in node.cigar.iter() {
                        match op {
                            Op::Match(l) => {
                                ops.extend(std::iter::repeat(kiley::Op::Match).take(*l))
                            }
                            Op::Del(l) => ops.extend(std::iter::repeat(kiley::Op::Del).take(*l)),
                            Op::Ins(l) => ops.extend(std::iter::repeat(kiley::Op::Ins).take(*l)),
                        }
                    }
                    (node.seq(), ops)
                })
                .unzip();
            let cov = seqs.len();
            let band_width = read_type.band_width(ref_unit.seq().len()) / 2;
            let start = std::time::Instant::now();
            let (consensus, asn, pss, score, k) = local_clustering_unit(
                ref_unit, &seqs, &mut ops, band_width, &mut rng, coverage, read_type, is_before,
            );
            for (node, ps) in units.iter_mut().zip(pss) {
                node.posterior = ps;
            }
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            // Change cigar string...
            let end = std::time::Instant::now();
            let elapsed = (end - start).as_millis();
            let len = consensus.len();
            println!(
                "ELAPSE\t{}\t{}\t{}\t{}\t{}\t{}",
                unit_id, elapsed, len, cov, k, is_before
            );
            (unit_id, (consensus, score, k))
        })
        .collect();
    for unit in ds.selected_chunks.iter_mut() {
        if let Some((consensus, score, cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = String::from_utf8(consensus.to_vec()).unwrap();
            unit.score = *score;
            unit.cluster_num = *cluster_num as usize;
        }
    }
}

type CLResult = (Vec<u8>, Vec<u8>, Vec<Vec<f64>>, f64, u8);
use kiley::Op;
fn local_clustering_unit<R: Rng>(
    ref_unit: &Unit,
    seqs: &[&[u8]],
    ops: &mut [Vec<Op>],
    band_width: usize,
    rng: &mut R,
    coverage: f64,
    read_type: definitions::ReadType,
    is_before: bool,
) -> CLResult {
    use haplotyper::local_clustering::kmeans;
    use haplotyper::local_clustering::kmeans::ClusteringConfig;
    let copy_num = ref_unit.copy_num as u8;
    let post = vec![vec![0f64]; seqs.len()];
    if is_before {
        let config = ClusteringConfig::new(band_width, copy_num, coverage, read_type);
        let hmm = kiley::gphmm::GPHMM::ont();
        let consensus =
            haplotyper::local_clustering::take_consensus(ref_unit, seqs, band_width, &hmm);
        if 1 < copy_num {
            let (asn, post, lk, k) =
                kmeans::clustering_dev(&consensus, &seqs, rng, &hmm, &config).unwrap();
            (consensus, asn, post, lk, k)
        } else {
            (consensus, vec![0; seqs.len()], post, 0f64, 1)
        }
    } else {
        let config = ClusteringConfig::new(band_width / 2, copy_num, coverage, read_type);
        let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
        let consensus = hmm.polish_until_converge_with(ref_unit.seq(), &seqs, ops, band_width);
        if 1 < ref_unit.copy_num {
            let (asn, post, lk, k) =
                kmeans::clustering_neo(&consensus, &seqs, ops, rng, &hmm, &config).unwrap();
            (consensus, asn, post, lk, k)
        } else {
            (consensus, vec![0; seqs.len()], post, 0f64, 1)
        }
    }
}

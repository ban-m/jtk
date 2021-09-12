#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum ExtractTarget {
    RawReads,
    HiCReads,
    Units,
    Assignments,
}

pub trait Extract {
    fn extract_fasta(&self, target: ExtractTarget) -> Vec<fasta::Record>;
    fn extract_assignments(&self) -> Vec<(usize, String, String)>;
}

use bio_utils::fasta;
impl Extract for definitions::DataSet {
    fn extract_assignments(&self) -> Vec<(usize, String, String)> {
        use std::collections::HashMap;
        let id2name: HashMap<_, _> = self
            .raw_reads
            .iter()
            .map(|r| (r.id, (r.name.clone(), r.desc.clone())))
            .collect();
        self.assignments
            .iter()
            .filter_map(|asn| {
                id2name
                    .get(&asn.id)
                    .map(|(n, d)| (asn.cluster, n.clone(), d.clone()))
            })
            .collect()
    }
    fn extract_fasta(&self, target: ExtractTarget) -> Vec<fasta::Record> {
        match target {
            ExtractTarget::RawReads => self
                .raw_reads
                .iter()
                .map(|r| {
                    let desc = if r.desc.is_empty() {
                        None
                    } else {
                        Some(r.desc.clone())
                    };
                    fasta::Record::with_data(&r.name, &desc, r.seq())
                })
                .collect(),
            ExtractTarget::HiCReads => self
                .hic_pairs
                .iter()
                .flat_map(|hic| {
                    let id1 = format!("{}_1", hic.pair_id);
                    let seq1 = fasta::Record::with_data(&id1, &None, hic.seq1());
                    let id2 = format!("{}_2", hic.pair_id);
                    let seq2 = fasta::Record::with_data(&id2, &None, hic.seq2());
                    vec![seq1, seq2]
                })
                .collect(),
            ExtractTarget::Units => self
                .selected_chunks
                .iter()
                .map(|u| {
                    let id = format!("{}", u.id);
                    fasta::Record::with_data(&id, &None, u.seq())
                })
                .collect(),
            _ => unreachable!(),
        }
    }
}

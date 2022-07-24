use definitions::*;
use kiley::hmm::guided::PairHiddenMarkovModel;
use std::collections::HashMap;
pub fn get_model(ds: &DataSet) -> Option<kiley::hmm::guided::PairHiddenMarkovModel> {
    ds.model_param.as_ref().map(|param| {
        let &HMMParam {
            mat_mat,
            mat_ins,
            mat_del,
            ins_mat,
            ins_ins,
            ins_del,
            del_mat,
            del_ins,
            del_del,
            ref mat_emit,
        } = param;
        PairHiddenMarkovModel {
            mat_mat,
            mat_ins,
            mat_del,
            ins_mat,
            ins_ins,
            ins_del,
            del_mat,
            del_ins,
            del_del,
            mat_emit: *mat_emit,
        }
    })
}

pub fn update_model(ds: &mut DataSet) {
    // let mut pileups: HashMap<u64, Vec<_>> =
    //     ds.selected_chunks.iter().map(|u| (u.id, vec![])).collect();
    let mut pileups: HashMap<_, Vec<_>> = HashMap::new();
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        pileups
            .entry((node.unit, node.cluster))
            .or_default()
            .push(node);
    }
    let hmm = estimate_model_parameters(ds.read_type, &pileups, &chunks);
    let PairHiddenMarkovModel {
        mat_mat,
        mat_ins,
        mat_del,
        ins_mat,
        ins_ins,
        ins_del,
        del_mat,
        del_ins,
        del_del,
        mat_emit,
    } = hmm;
    ds.model_param = Some(HMMParam {
        mat_mat,
        mat_ins,
        mat_del,
        ins_mat,
        ins_ins,
        ins_del,
        del_mat,
        del_ins,
        del_del,
        mat_emit,
    });
}

fn estimate_model_parameters<N: std::borrow::Borrow<Node>>(
    read_type: ReadType,
    pileups: &HashMap<(u64, u64), Vec<N>>,
    chunks: &HashMap<u64, &Unit>,
) -> kiley::hmm::guided::PairHiddenMarkovModel {
    let mut covs: Vec<_> = pileups.iter().map(|x| x.1.len()).collect();
    let (_, &mut cov, _) = covs.select_nth_unstable(pileups.len() / 2);
    let mut seqs_and_ref_units: Vec<_> = pileups
        .iter()
        .filter(|(_, us)| (cov.max(2) - 2..cov + 2).contains(&us.len()))
        .map(|((unit, _), us)| (chunks.get(unit).unwrap(), us))
        .collect();
    seqs_and_ref_units.sort_by_cached_key(|c| c.0.id);
    seqs_and_ref_units.truncate(2);
    for (chunk, units) in seqs_and_ref_units.iter() {
        debug!("LOCAL\tSAMPLE\t{}\t{}", chunk.id, units.len());
    }
    let mut hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let mut polishing_pairs: Vec<_> = seqs_and_ref_units
        .iter()
        .map(|(ref_unit, nodes)| {
            let band_width = read_type.band_width(ref_unit.seq().len());
            let ops: Vec<Vec<_>> = nodes
                .iter()
                .map(|n| crate::misc::ops_to_kiley(&n.borrow().cigar))
                .collect();
            let seqs: Vec<_> = nodes.iter().map(|n| n.borrow().seq()).collect();
            (ref_unit.seq().to_vec(), seqs, ops, band_width)
        })
        .collect();
    // polishing_pairs
    //     .par_iter_mut()
    //     .for_each(|(consensus, seqs, ops, bw)| {
    //         use kiley::bialignment::guided;
    //         *consensus = guided::polish_until_converge_with(consensus, seqs, ops, *bw);
    //         *consensus = hmm.polish_until_converge_with(consensus, seqs, ops, *bw);
    //     });
    debug!("TUNING");
    for _ in 0..3 {
        for (consensus, seqs, ops, bw) in polishing_pairs.iter_mut() {
            use kiley::bialignment::guided;
            *consensus = guided::polish_until_converge_with(consensus, seqs, ops, *bw);
            *consensus = hmm.polish_until_converge_with(consensus, seqs, ops, *bw);
            hmm.fit_naive_with_par(consensus, seqs, ops, *bw);
        }
    }
    debug!("HMM\n{}", hmm);
    hmm
}

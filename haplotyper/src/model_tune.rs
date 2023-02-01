use definitions::*;
use kiley::hmm::guided::HMMConfig;
use std::collections::HashMap;

use kiley::hmm::PairHiddenMarkovModel;
use kiley::hmm::PairHiddenMarkovModelOnStrands;

pub trait ModelFit {
    fn fit_models_on_both_strands(&self) -> Option<PairHiddenMarkovModelOnStrands>;
    fn get_model(&self) -> Option<PairHiddenMarkovModel>;
}

impl ModelFit for DataSet {
    fn fit_models_on_both_strands(&self) -> Option<PairHiddenMarkovModelOnStrands> {
        estimate_model_parameters_on_both_strands(self)
    }
    fn get_model(&self) -> Option<PairHiddenMarkovModel> {
        get_model(self)
    }
}

pub fn get_model(ds: &DataSet) -> Option<kiley::hmm::PairHiddenMarkovModel> {
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
            ref ins_emit,
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
            ins_emit: *ins_emit,
        }
    })
}

pub fn update_model(ds: &mut DataSet) {
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
        ins_emit,
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
        ins_emit,
    });
}

const TRAIN_UNIT_SIZE: usize = 3;
const TRAIN_ROUND: usize = 3;
fn estimate_model_parameters_on_both_strands(
    ds: &DataSet,
) -> Option<PairHiddenMarkovModelOnStrands> {
    fn truncate_bucket<T>(pileups: HashMap<u64, Vec<T>>) -> Vec<(u64, Vec<T>)> {
        let mut covs: Vec<_> = pileups.values().map(|x| x.len()).collect();
        let (_, &mut cov, _) = covs.select_nth_unstable(pileups.len() / 2);
        pileups
            .into_iter()
            .filter(|(_, us)| (cov.max(2) - 2..cov + 2).contains(&us.len()))
            .take(TRAIN_UNIT_SIZE)
            .collect()
    }
    let model = ds.get_model()?;
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let models = PairHiddenMarkovModelOnStrands::new(model.clone(), model);
    let pileups = {
        let mut pileups: HashMap<_, _> = chunks.keys().map(|&k| (k, vec![])).collect();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            pileups.entry(node.unit).or_default().push(node);
        }
        truncate_bucket(pileups)
    };
    let models = train(ds.read_type, &chunks, &pileups, &models);
    fn train(
        read_type: ReadType,
        chunks: &HashMap<u64, &Unit>,
        pileups: &[(u64, Vec<&Node>)],
        models: &PairHiddenMarkovModelOnStrands,
    ) -> PairHiddenMarkovModelOnStrands {
        let mut models = models.clone();
        let mut polishing_pairs: Vec<_> = pileups
            .iter()
            .filter_map(|(uid, nodes)| {
                let ref_unit = chunks.get(uid)?;
                let band_width = read_type.band_width(ref_unit.seq().len());
                let ops: Vec<Vec<_>> = nodes
                    .iter()
                    .map(|n| crate::misc::ops_to_kiley(&n.cigar))
                    .collect();
                let seqs: Vec<_> = nodes.iter().map(|n| n.seq()).collect();
                let strands: Vec<_> = nodes.iter().map(|n| n.is_forward).collect();
                let chunk_seq = ref_unit.seq().to_vec();
                Some((chunk_seq, seqs, ops, strands, band_width))
            })
            .collect();
        assert!(!polishing_pairs.is_empty());
        use rayon::prelude::*;
        let bw = polishing_pairs.iter().map(|x| x.4).max().unwrap();
        for _ in 0..TRAIN_ROUND {
            polishing_pairs
                .par_iter_mut()
                .for_each(|(cons, seqs, ops, strands, bw)| {
                    let config = HMMConfig::new(*bw, seqs.len(), 0);
                    *cons =
                        models.polish_until_converge_with_conf(cons, seqs, ops, strands, &config);
                });
            let training_datapack: Vec<_> = polishing_pairs
                .iter()
                .map(|(cons, seqs, ops, strands, _)| {
                    kiley::hmm::TrainingDataPack::new(cons, strands, seqs, ops)
                })
                .collect();
            models.fit_multiple_with_par(&training_datapack, bw);
        }
        debug!("FORWARD\tHMM\n{}", models.forward());
        debug!("REVERSE\tHMM\n{}", models.reverse());
        models.clone()
    }
    Some(models)
}

fn estimate_model_parameters<N: std::borrow::Borrow<Node>>(
    read_type: ReadType,
    pileups: &HashMap<(u64, u64), Vec<N>>,
    chunks: &HashMap<u64, &Unit>,
) -> kiley::hmm::PairHiddenMarkovModel {
    let mut covs: Vec<_> = pileups.iter().map(|x| x.1.len()).collect();
    let (_, &mut cov, _) = covs.select_nth_unstable(pileups.len() / 2);
    let mut seqs_and_ref_units: Vec<_> = pileups
        .iter()
        .filter(|(_, us)| (cov.max(2) - 2..cov + 2).contains(&us.len()))
        .map(|((unit, _), us)| (chunks.get(unit).unwrap(), us))
        .collect();
    seqs_and_ref_units.sort_by_cached_key(|c| c.0.id);
    seqs_and_ref_units.truncate(TRAIN_UNIT_SIZE);
    for (chunk, units) in seqs_and_ref_units.iter() {
        debug!("LOCAL\tSAMPLE\t{}\t{}", chunk.id, units.len());
    }
    let mut hmm = kiley::hmm::PairHiddenMarkovModel::default();
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
    for _ in 0..TRAIN_ROUND {
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

// use kiley::Op;
// const PRIOR: f64 = 10f64;
// pub fn fine_tune(
//     hmm: &PairHiddenMarkovModel,
//     template: &[u8],
//     seqs: &[&[u8]],
//     ops: &[Vec<Op>],
// ) -> PairHiddenMarkovModel {
//     let PairHiddenMarkovModel {
//         mat_mat,
//         mat_ins,
//         mat_del,
//         ins_mat,
//         ins_ins,
//         ins_del,
//         del_mat,
//         del_ins,
//         del_del,
//         mut mat_emit,
//         mut ins_emit,
//     } = hmm.clone();
//     let mut transitions = [
//         [mat_mat, mat_ins, mat_del],
//         [ins_mat, ins_ins, ins_del],
//         [del_mat, del_ins, del_del],
//     ];
//     transitions.iter_mut().flatten().for_each(|x| *x *= PRIOR);
//     mat_emit.iter_mut().for_each(|x| *x *= PRIOR);
//     ins_emit.iter_mut().for_each(|x| *x *= PRIOR);
//     for (seq, ops) in seqs.iter().zip(ops.iter()) {
//         let params = (&mut transitions, &mut mat_emit, &mut ins_emit);
//         register_alignments(template, seq, ops, params);
//     }
//     let mat = (transitions[0][0], transitions[0][1], transitions[0][2]);
//     let ins = (transitions[1][0], transitions[1][1], transitions[1][2]);
//     let del = (transitions[2][0], transitions[2][1], transitions[2][2]);
//     PairHiddenMarkovModel::new(mat, ins, del, &mat_emit, &ins_emit)
// }

// pub fn register_alignments(
//     template: &[u8],
//     xs: &[u8],
//     ops: &[Op],
//     (transitions, mat_emit, ins_emit): (&mut [[f64; 3]; 3], &mut [f64; 16], &mut [f64; 20]),
// ) {
//     fn op_to_state(op: Op) -> usize {
//         match op {
//             Op::Mismatch => 0,
//             Op::Match => 0,
//             Op::Ins => 1,
//             Op::Del => 2,
//         }
//     }
//     fn base_to_idx(base: u8) -> usize {
//         match base {
//             b'A' | b'a' => 0,
//             b'C' | b'c' => 1,
//             b'G' | b'g' => 2,
//             b'T' | b't' => 3,
//             _ => panic!(),
//         }
//     }
//     let mut state = op_to_state(ops[0]);
//     let (mut rpos, mut qpos) = (0, 0);
//     let rbase = base_to_idx(template[rpos]) << 2;
//     let qbase = base_to_idx(xs[qpos]);
//     match state {
//         0 => {
//             mat_emit[rbase | qbase] += 1f64;
//             rpos += 1;
//             qpos += 1;
//         }
//         1 => {
//             ins_emit[rbase | qbase] += 1f64;
//             qpos += 1;
//         }
//         _ => {
//             rpos += 1;
//         }
//     }
//     for op in ops.iter().skip(1) {
//         let next = op_to_state(*op);
//         transitions[state][next] += 1f64;
//         state = next;
//         if state == 2 {
//             rpos += 1;
//             continue;
//         }
//         let rbase = base_to_idx(template[rpos]) << 2;
//         let qbase = base_to_idx(xs[qpos]);
//         match state {
//             0 => {
//                 mat_emit[rbase | qbase] += 1f64;
//                 rpos += 1;
//                 qpos += 1;
//             }
//             _ => {
//                 ins_emit[rbase | qbase] += 1f64;
//                 qpos += 1;
//             }
//         }
//     }
//     for &base in template.iter() {
//         ins_emit[16 + base_to_idx(base)] += 1.0;
//     }
// }

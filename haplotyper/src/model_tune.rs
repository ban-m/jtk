use definitions::*;
use kiley::hmm::guided::HMMConfig;
use std::collections::HashMap;

use kiley::hmm::PairHiddenMarkovModel;
use kiley::hmm::PairHiddenMarkovModelOnStrands;

pub trait ModelFit {
    fn fit_models_on_both_strands(&self) -> Option<PairHiddenMarkovModelOnStrands>;
    fn update_models_on_both_strands(&mut self);
    fn get_model(&self) -> PairHiddenMarkovModel;
}

impl ModelFit for DataSet {
    fn fit_models_on_both_strands(&self) -> Option<PairHiddenMarkovModelOnStrands> {
        estimate_model_parameters_on_both_strands(self)
    }
    fn update_models_on_both_strands(&mut self) {
        let model = self.fit_models_on_both_strands().unwrap();
        let forward = kiley_into_def(model.forward());
        let reverse = kiley_into_def(model.reverse());
        self.model_param = HMMParamOnStrands { forward, reverse };
    }
    fn get_model(&self) -> PairHiddenMarkovModel {
        get_model(self)
    }
}

fn kiley_into_def(model: &kiley::hmm::PairHiddenMarkovModel) -> definitions::HMMParam {
    let &PairHiddenMarkovModel {
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
    } = model;
    HMMParam {
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
}

fn get_model(ds: &DataSet) -> kiley::hmm::PairHiddenMarkovModel {
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
    } = &ds.model_param.forward;
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
    let model = ds.get_model();
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
        chunks: &HashMap<u64, &Chunk>,
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

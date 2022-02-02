//! Tiny implementation of a string graph.
//! Implementation remarks:
//! - It does/should not implement consensus generation procedure so seriously. It should
//!   output runs of (unit,cluster)s.
//! - It is too much for a string graph to rememeber all the reads, including contained ones.
//!   In fact, if we allow the graph to 'forget' about reads contained by others, it would be
//!   much easier to implement.
use definitions::*;
use rayon::prelude::*;

use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct AssembleConfig {
    identity: f64,
    ovlp: i32,
}

impl AssembleConfig {
    pub fn new(ovlp: i32, identity: f64) -> Self {
        Self { identity, ovlp }
    }
}

const IDEN: f64 = 0.3;
const OVLP: i32 = 3;
pub const DEFAULT_CONFIG: AssembleConfig = AssembleConfig {
    identity: IDEN,
    ovlp: OVLP,
};

/// This trait requires some functions to "align" reads to each other.
/// It is OK to return None to indicate that two sequence can not be aligned to each other.
/// Also, some informations can be embeded by using alignment configuration trait.
pub trait Alignable: Sized {
    type Config: AlignmentConfigure;
    fn align(x1: &[Self], x2: &[Self], conf: &Self::Config) -> Option<Alignment>;
}

pub trait AlignmentConfigure {}

#[derive(Clone)]
pub struct Alignment {
    // The correct probability. shoule be in 0..1f64.
    pub correct_prob: f64,
    pub ovlp_len_from: usize,
    pub ovlp_len_to: usize,
    pub alntype: AlnType,
}

impl Alignment {
    pub fn new(
        correct_prob: f64,
        ovlp_len_from: usize,
        ovlp_len_to: usize,
        alntype: AlnType,
    ) -> Self {
        Self {
            correct_prob,
            ovlp_len_from,
            ovlp_len_to,
            alntype,
        }
    }
}

#[derive(Clone)]
pub enum AlnType {
    /// -------->
    ///     --------->
    TailToHead,
    /// --------->
    ///     <-----------
    TailToTail,
    /// <-----------
    ///      ------------>
    HeadToHead,
    /// <------------
    ///        <---------
    HeadToTail,
    /// ------------->
    ///       ---->
    Containing,
    ///      ----->
    /// -------------->
    Contained,
    ///    ------->
    ///    ------->
    Diagonal,
}

#[derive(Clone)]
pub struct StringGraph<T: Clone + Alignable> {
    // nodes of string graph.
    pub nodes: Vec<Option<SGNode<T>>>,
}

#[derive(Clone)]
pub struct SGNode<T: Clone + Alignable> {
    edges: Vec<SGEdge>,
    // Sequence inside it.
    seq: Vec<T>,
}

#[derive(Clone)]
pub struct SGEdge {
    from: usize,
    to: usize,
    alignment: Alignment,
}

fn match_score(u1: (u64, bool), u2: (u64, bool)) -> i32 {
    if u1 == u2 {
        1
    } else {
        -1000
    }
}

#[derive(Debug, Clone)]
struct NodeAlignConfig {
    prior: f64,
}

impl AlignmentConfigure for NodeAlignConfig {}

type DirNode<'a> = ((u64, bool), &'a [f64]);
impl<'a> Alignable for DirNode<'a> {
    type Config = NodeAlignConfig;
    fn align(x1: &[Self], x2: &[Self], c: &Self::Config) -> Option<Alignment> {
        let forward = align_forward(x1, x2, c);
        let backward = align_reverse(x1, x2, c);
        match (forward, backward) {
            (Some(f), Some(b)) if f.0 < b.0 => todo!(),
            (Some(f), Some(b)) => todo!(),
            (Some(f), None) => todo!(),
            (None, Some(b)) => todo!(),
            _ => None,
        }
    }
}

fn align_reverse<'a>(
    x1: &[DirNode<'a>],
    x2: &[DirNode<'a>],
    c: &NodeAlignConfig,
) -> Option<(f64, usize, usize, usize, usize)> {
    let mut x2 = x2.to_vec();
    x2.reverse();
    x2.iter_mut().for_each(|x| x.0 .1 = !x.0 .1);
    let len = x2.len();
    align_forward(x1, &x2, c).map(|(prob, s1, e1, s2, e2)| (prob, s1, e1, len - e2, len - s2))
}
fn align_forward<'a>(
    x1: &[DirNode<'a>],
    x2: &[DirNode<'a>],
    c: &NodeAlignConfig,
) -> Option<(f64, usize, usize, usize, usize)> {
    // First, compute the alignment path,
    // which is almost unique.
    let mut dp = vec![vec![0; x2.len() + 1]; x1.len() + 1];
    for (i, (u1, _)) in x1.iter().enumerate() {
        let i = i + 1;
        for (j, (u2, _)) in x2.iter().enumerate() {
            let j = j + 1;
            let mat = match_score(*u1, *u2);
            dp[i][j] = (dp[i][j - 1] - 1).max(dp[i - 1][j] - 1).max(dp[i][j] + mat);
        }
    }
    let rowmax = dp.last().unwrap().iter().enumerate().max_by_key(|x| x.1);
    let (rowidx, rowmax) = rowmax.unwrap();
    let colmax = dp.iter().filter_map(|x| x.last()).enumerate();
    let (colidx, colmax) = colmax.max_by_key(|x| x.1).unwrap();
    let (x1end, x2end) = if rowmax < colmax {
        (colidx, x2.len())
    } else {
        (x1.len(), rowidx)
    };
    if x1end == 0 || x2end == 0 {
        return None;
    }
    // Construct alignment path.
    let (mut i, mut j) = (x1end, x2end);
    let mut ops = vec![];
    while 0 < i || 0 < j {
        let score = dp[i][j];
        if score == dp[i][j - 1] - 1 {
            ops.push(2);
            j -= 1;
        } else if score == dp[i - 1][j] - 1 {
            ops.push(1);
            i -= 1;
        } else {
            ops.push(0);
            let mat = match_score(x1[i - 1].0, x2[j - 1].0);
            assert_eq!(score, dp[i - 1][j - 1] + mat);
            i -= 1;
            j -= 1;
        }
    }
    let (x1start, x2start) = (i, j);
    ops.reverse();
    let prob = get_posterior(&x1[x1start..x1end], &x2[x2start..x2end], &ops);
    Some((prob, x1start, x1end, x2start, x2end))
}

fn get_posterior<'a>(x1: &[DirNode<'a>], x2: &[DirNode<'a>], ops: &[u8]) -> f64 {
    todo!()
}

impl<'a> StringGraph<DirNode<'a>> {
    //! Construct a new string graph.
    pub fn new(dataset: &DataSet, config: &AssembleConfig) -> Self {
        todo!()
    }
    pub fn clean_up(&mut self, config: &AssembleConfig) {
        todo!()
    }
    pub fn spell(&self) -> Vec<Vec<(u64, u64, bool)>> {
        todo!()
    }
}

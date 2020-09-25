pub fn consensus(seq: &[Vec<u8>]) -> Vec<u8> {
    if seq.is_empty() {
        return vec![];
    }
    let mut seqs: Vec<_> = seq.iter().map(|e| e.as_slice()).collect();
    seqs.sort_by_key(|x| x.len());
    let template = seqs.pop().unwrap();
    seqs.iter()
        .fold(poa_hmm::POA::new(template, 1.), |x, y| {
            x.add_default_banded(y, 150)
        })
        .consensus()
    // let alignments: Vec<_> = seqs
    //     .par_iter()
    //     .map(|q| edlib_sys::local(template, q))
    //     .collect();
    // let pileup = seqs
    //     .iter()
    //     .zip(alignments.iter())
    //     .fold(PileUp::new(template), |x, (seq, y)| x.add(seq, y));
    // let mut result = vec![];
    // for column in pileup.columns.iter() {
    //     column.generate(&mut result);
    // }
    // result
}

// #[derive(Debug, Clone)]
// struct PileUp {
//     columns: Vec<Column>,
// }

// impl PileUp {
//     fn new(seq: &[u8]) -> Self {
//         let columns: Vec<_> = seq.iter().map(|&b| Column::new(b)).collect();
//         Self { columns }
//     }
// }

// #[derive(Debug, Clone)]
// struct Column {
//     coverage: u32,
//     matches: [u8; 4],
//     insertions: [u8; 4],
//     deletions: u32,
// }

// fn b2i(b: u8) -> usize {
//     match b {
//         b'A' => 0,
//         b'C' => 1,
//         b'G' => 2,
//         b'T' => 3,
//         _ => panic!(),
//     }
// }

// fn i2b(i: usize) -> u8 {
//     match i {
//         0 => b'A',
//         1 => b'C',
//         2 => b'G',
//         3 => b'T',
//         _ => panic!(),
//     }
// }

// impl Column {
//     fn new(b: u8) -> Self {
//         let mut matches = [0; 4];
//         matches[b2i(b)] += 1;
//         Self {
//             coverage: 1,
//             matches,
//             insertions: [0; 4],
//             deletions: 0,
//         }
//     }
//     fn generate(&self, buf:&mut Vec<u8>){

//     }
// }

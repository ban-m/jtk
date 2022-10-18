fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let references = bio_utils::fasta::parse_into_vec(&args[1])?;
    let assemblies = bio_utils::fasta::parse_into_vec(&args[2])?;
    let ref_hap1: Vec<_> = args[3].split_whitespace().collect();
    let ref_hap2: Vec<_> = args[4].split_whitespace().collect();
    let asm_hap1: Vec<_> = args[5].split_whitespace().collect();
    let asm_hap2: Vec<_> = args[6].split_whitespace().collect();
    let filter = args
        .get(7)
        .map(|x| x.parse::<usize>().unwrap())
        .unwrap_or(0);
    let sum: usize = assemblies
        .iter()
        .filter(|asm| {
            asm_hap1
                .iter()
                .chain(asm_hap2.iter())
                .any(|&name| name == asm.id())
        })
        .map(|record| record.seq().len())
        .sum();
    let num = asm_hap1.len() + asm_hap2.len();
    // cat "$ASM" | paste - - | awk '(length($2) > 500000){sum+=1;total+=length($2)}END{printf(" & %d & %d & ",sum,int(total/1000))}'
    let (dist1, cov1) = {
        // h1 - h1, h2 - h2.
        let (dist_h11, cov11) = compare_dev(&references, &ref_hap1, &assemblies, &asm_hap1, filter);
        let (dist_h22, cov22) = compare_dev(&references, &ref_hap2, &assemblies, &asm_hap2, filter);
        (dist_h11 + dist_h22, cov11 + cov22)
    };
    eprintln!("{dist1}\t{cov1}");

    let (dist2, cov2) = {
        let (dist_h12, cov12) = compare_dev(&references, &ref_hap1, &assemblies, &asm_hap2, filter);
        let (dist_h21, cov21) = compare_dev(&references, &ref_hap2, &assemblies, &asm_hap1, filter);
        (dist_h12 + dist_h21, cov12 + cov21)
    };
    eprintln!("{dist2}\t{cov2}");
    // let dist = dist1.min(dist2);
    // Converts to QV ...
    let error_rate = if dist1 < dist2 {
        dist1 as f64 / cov1 as f64
    } else {
        dist2 as f64 / cov2 as f64
    };
    let qv = -10f64 * error_rate.log10();
    println!(" & {num} & {} & {qv:.1}", sum / 1000);
    Ok(())
}

fn compare_dev(
    refers: &[Record],
    ref_hap: &[&str],
    asms: &[Record],
    asm_hap: &[&str],
    filter: usize,
) -> (usize, usize) {
    use rand::thread_rng;
    use rand::Rng;
    use std::io::Write;
    let mut rng = thread_rng();
    let q_file_name = format!("{}.fa", rng.gen_range(0..1000));
    let r_file_name = format!("{}.fa", rng.gen_range(0..1000));
    {
        let mut q_file = std::fs::File::create(&q_file_name)
            .map(std::io::BufWriter::new)
            .unwrap();
        let mut r_file = std::fs::File::create(&r_file_name)
            .map(std::io::BufWriter::new)
            .unwrap();
        for name in asm_hap {
            let asm = asms.iter().find(|r| r.id() == *name).unwrap();
            writeln!(&mut q_file, "{asm}").unwrap();
        }
        for ref_name in ref_hap.iter() {
            let record = refers.iter().find(|r| r.id() == *ref_name).unwrap();
            writeln!(&mut r_file, "{record}").unwrap();
        }
    }
    let arg = ["-c", "--secondary=no", "--eqx", "-x", "asm20"];
    let output = haplotyper::minimap2::minimap2_args(&q_file_name, &r_file_name, &arg);
    std::fs::remove_file(&q_file_name).unwrap();
    std::fs::remove_file(&r_file_name).unwrap();
    String::from_utf8_lossy(&output)
        .lines()
        .map(|line| {
            for field in line.split_whitespace() {
                if field.len() < 50 {
                    eprint!("{field}\t");
                }
            }
            eprintln!();
            let paf = bio_utils::paf::PAF::new(line).unwrap();
            let dist = paf.get_tag("NM").unwrap().1.parse::<usize>().unwrap();
            let covered = paf.blocklen;
            if filter < covered {
                (dist, covered)
            } else {
                (0, 0)
            }
        })
        .fold((0, 0), |(d, cov), (x, y)| (d + x, cov + y))
}

use bio_utils::fasta::Record;
// fn compare(refers: &[Record], ref_hap: &[&str], asms: &[Record], asm_hap: &[&str]) -> usize {
//     let asms: Vec<_> = asm_hap
//         .iter()
//         .filter_map(|name| asms.iter().find(|r| r.id() == *name))
//         .collect();
//     ref_hap
//         .iter()
//         .map(|ref_name| {
//             let record = refers.iter().find(|r| r.id() == *ref_name).unwrap();
//             asms.iter()
//                 .map(|asm| {
//                     let (query, target) = match asm.seq().len() < record.seq().len() {
//                         true => (asm.seq(), record.seq()),
//                         false => (record.seq(), asm.seq()),
//                     };
//                     if similar(query, target) {
//                         dist(query, target)
//                     } else {
//                         let query = bio_utils::revcmp(query);
//                         if similar(&query, target) {
//                             dist(&query, target)
//                         } else {
//                             query.len() + target.len()
//                         }
//                     }
//                 })
//                 .min()
//                 .unwrap()
//         })
//         .sum()
// }

// fn dist(query: &[u8], refr: &[u8]) -> usize {
//     let mode = edlib_sys::AlignMode::Global;
//     let task = edlib_sys::AlignTask::Alignment;
//     let aln = edlib_sys::align(&query, &refr, mode, task);
//     let start = std::time::Instant::now();
//     let path = haplotyper::misc::edlib_to_kiley(&aln.operations().unwrap());
//     let edl = std::time::Instant::now();
//     let path =
//         kiley::bialignment::guided::global_guided(&refr, &query, &path, 5, (1, -1, -1, -1)).1;
//     let path =
//         kiley::bialignment::guided::global_guided(&refr, &query, &path, 5, (1, -1, -1, -1)).1;
//     let path = remove_indel(path);
//     let fix = std::time::Instant::now();
//     let edlib = (edl - start).as_secs();
//     let fixing = (fix - edl).as_secs();
//     eprintln!("{edlib}\t{fixing}");
//     path.iter().filter(|&&op| op != kiley::Op::Match).count()
// }

// fn remove_indel(mut path: Vec<kiley::Op>) -> Vec<kiley::Op> {
//     while path.last() == Some(&kiley::Op::Del) || path.last() == Some(&kiley::Op::Ins) {
//         path.pop();
//     }
//     path.reverse();
//     while path.last() == Some(&kiley::Op::Del) || path.last() == Some(&kiley::Op::Ins) {
//         path.pop();
//     }
//     path.reverse();
//     path
// }

// fn similar(xs: &[u8], ys: &[u8]) -> bool {
//     use std::collections::HashSet;
//     let xs: HashSet<_> = xs.windows(15).map(|x| x.to_vec()).collect();
//     let ys: HashSet<_> = ys.windows(15).map(|x| x.to_vec()).collect();
//     let intersection = xs.intersection(&ys).count();
//     let union = xs.union(&ys).count();
//     eprintln!("{intersection},{union}");
//     (intersection as f64 / union as f64) > 0.05
// }

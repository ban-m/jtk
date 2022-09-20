fn main() -> std::io::Result<()> {
    env_logger::init();
    let mode = edlib_sys::AlignMode::Global;
    let task = edlib_sys::AlignTask::Alignment;
    for len in [100, 1000, 10_000, 100_000, 1_000_000] {
        let mut rng = rand::thread_rng();
        let seq = kiley::gen_seq::generate_seq(&mut rng, len);
        let prof = kiley::gen_seq::Profile::new(0.001, 0.001, 0.001);
        let seq2 = kiley::gen_seq::introduce_randomness(&seq, &mut rng, &prof);
        let start = std::time::Instant::now();
        let aln = edlib_sys::align(&seq, &seq2, mode, task);
        let path = haplotyper::misc::edlib_to_kiley(&aln.operations().unwrap());
        let path =
            kiley::bialignment::guided::global_guided(&seq2, &seq, &path, 10, (1, -1, -1, -1)).1;
        let dist = path.iter().filter(|&&op| op != kiley::Op::Match).count();
        let end = std::time::Instant::now();
        let time = (end - start).as_millis();
        println!("{dist}\t{time}");
    }
    Ok(())
}

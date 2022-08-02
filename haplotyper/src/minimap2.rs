/// Invole minimap2 with specified input and return the stdout directry.
pub fn minimap2(
    query: &str,
    target: &str,
    threads: usize,
    preset: &str,
    is_ava: bool,
    is_sam: bool,
) -> Vec<u8> {
    let thr = format!("{}", threads);
    let mut args = vec!["-x", preset];
    if is_ava {
        args.push("-X");
    }
    if is_sam {
        args.push("-a")
    } else {
        args.push("-c")
    }
    args.extend(&["-t", &thr, target, query]);
    let aln = std::process::Command::new("minimap2")
        .args(&args)
        .output()
        .unwrap();
    if !aln.status.success() {
        panic!("Minimap2,{:?}", String::from_utf8_lossy(&aln.stderr));
    } else {
        aln.stdout
    }
}

pub fn minimap2_args(target: &str, query: &str, args: &[&str]) -> Vec<u8> {
    let mut args = args.to_vec();
    args.push(target);
    args.push(query);
    let aln = std::process::Command::new("minimap2").args(&args).output();
    let aln = match aln {
        Ok(res) => res,
        Err(why) => {
            eprintln!("{target},{query},{args:?}");
            eprintln!("{why:?}");
            panic!()
        }
    };
    if !aln.status.success() {
        panic!("Minimap2,{}", std::str::from_utf8(&aln.stderr).unwrap());
    } else {
        aln.stdout
    }
}

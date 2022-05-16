// Synopsis
// `cargo run --release --bin simplify_gfa_from_flye -- ${GFA from Flye} <SpellPattern>+ > output.fa`
// Where,
//    - <SpellPattern> = <ContigSpec>(:<ContiSpec>)*" "
//    - <ContigSpec> = ${NAME},{+/-}
// The ID(or the string ater >) would be automatically 0,1,2,3,....
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::io::*;
    let nodes: HashMap<_, _> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| l.starts_with("S\t"))
        .map(|line| {
            let mut line = line.split_whitespace();
            let id = line.nth(1).unwrap();
            let seq = line.next().unwrap();
            (id.to_string(), seq.to_string())
        })
        .collect();
    for (i, spec) in args[2..].iter().enumerate() {
        let seq: Vec<_> = spec
            .split(':')
            .flat_map(|ctgspec| {
                let (ctgname, dir) = ctgspec.split_once(',').unwrap();
                nodes
                    .get(ctgname)
                    .map(|seq| match dir {
                        "+" => seq.as_bytes().to_vec(),
                        "-" => bio_utils::revcmp(seq.as_bytes()),
                        _ => panic!("{}", dir),
                    })
                    .unwrap()
            })
            .collect();
        println!(">{i}");
        println!("{}", std::str::from_utf8(&seq).unwrap());
    }
    Ok(())
}

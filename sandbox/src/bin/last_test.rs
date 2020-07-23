fn main() {
    let args: Vec<_> = std::env::args().collect();
    let reference = &args[1];
    let reads = &args[2];
    let lastdb = std::process::Command::new("lastdb")
        .args(&["-R", "00", "-Q", "0", "tempref", reference])
        .status()
        .unwrap()
        .success();
    if !lastdb {
        panic!("something went wrong.")
    }
    let last_train = std::process::Command::new("last-train")
        .args(&["-P", "23", "-Q", "0", "tempref", reads])
        .output()
        .unwrap();
    use std::io::{BufWriter, Write};
    {
        let mut wtr = BufWriter::new(std::fs::File::create("param.par").unwrap());
        wtr.write_all(&last_train.stdout).unwrap();
        wtr.flush().unwrap();
    }
    let lastal = std::process::Command::new("lastal")
        .args(&[
            "-f",
            "tab",
            "-P",
            "23",
            "-R",
            "00",
            "-Q",
            "0",
            "-p",
            "param.par",
            "tempref",
            reads,
        ])
        .output()
        .unwrap();
    eprintln!("{}", lastal.status);
    eprintln!("{}", String::from_utf8_lossy(&lastal.stderr));
    println!("{}", String::from_utf8_lossy(&lastal.stdout));
}

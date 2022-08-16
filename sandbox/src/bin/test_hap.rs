fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    use std::io::*;
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use std::collections::HashSet;
    let selection: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    use haplotyper::local_clustering::*;
    local_clustering_selected(&mut ds, &selection);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}

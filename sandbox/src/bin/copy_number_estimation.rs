use std::io::BufReader;
fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::copy_number_estimation::*;
    let config = Config::default();
    let (node_cp, edge_cp) = ds.estimate_copy_numbers(&config);
    println!("NODE\tunit\tcluster\tcopy.number");
    for ((unit, cluster), cp) in node_cp {
        println!("NODE\t{}\t{}\t{}", unit, cluster, cp);
    }
    println!("EDGE\tfrom.unit\tfrom.cluster\tto.unit\tto.cluster\tcopy.number");
    for (((fu, fc), (tu, tc)), cp) in edge_cp {
        println!("EDGE\t{}\t{}\t{}\t{}\t{}", fu, fc, tu, tc, cp);
    }
}

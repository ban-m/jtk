use std::fs::File;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let file = File::open(&args[1])?;
    let gfa = gfa::GFA::from_reader(file);
    for record in gfa {
        println!("{}", record);
    }
    Ok(())
}

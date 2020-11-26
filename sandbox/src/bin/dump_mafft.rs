use std::fs::File;
use std::io::{BufReader, Read};
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut buf = String::new();
    let mut rdr = BufReader::new(File::open(&args[1]).unwrap());
    rdr.read_to_string(&mut buf).unwrap();
    for record in buf.split('>') {
        for line in record.split('\n').skip(1) {
            print!("{}", line);
        }
        println!();
    }
}

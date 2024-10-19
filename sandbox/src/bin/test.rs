use serde::Deserialize;
use serde::Serialize;
use toml::from_str;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Test {
    a_or_b: std::path::PathBuf,
}
fn main() -> std::io::Result<()> {
    let content = r#"
    a_or_b = "./test"
    "#;

    let file: Test = from_str(content).unwrap();
    println!("{file:?}");
    Ok(())
}

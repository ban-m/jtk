//! Definition of the position of the contig.
/// Position::Head is the up-stream position of a unit, and Tail is downstream by default.
#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash, PartialOrd, Ord)]
pub enum Position {
    Head,
    Tail,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = match self {
            Position::Head => 'H',
            Position::Tail => 'T',
        };
        write!(f, "{}", x)
    }
}

impl std::ops::Not for Position {
    type Output = Self;
    fn not(self) -> Self {
        match self {
            Position::Head => Position::Tail,
            Position::Tail => Position::Head,
        }
    }
}

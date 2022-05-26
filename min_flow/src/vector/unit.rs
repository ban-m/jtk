//!
//! Trait for types which have a unit of addition/multiplication
//!
use crate::prob::Prob;

///
/// get a unit of addition
///
pub trait UnitAdd: Copy + PartialEq {
    fn unit_add() -> Self;
    fn is_unit_add(&self) -> bool {
        *self == Self::unit_add()
    }
}

impl UnitAdd for u32 {
    fn unit_add() -> u32 {
        0
    }
}
impl UnitAdd for usize {
    fn unit_add() -> usize {
        0
    }
}
impl UnitAdd for f64 {
    fn unit_add() -> f64 {
        0.0
    }
}
impl UnitAdd for Prob {
    fn unit_add() -> Prob {
        Prob::from_log_prob(f64::NEG_INFINITY)
    }
}

///
/// get a unit of multiplication
///
pub trait UnitMul: Copy + PartialEq {
    fn unit_mul() -> Self;
    fn zero_mul() -> Self;
    fn is_unit_mul(&self) -> bool {
        *self == Self::unit_mul()
    }
    fn is_zero_mul(&self) -> bool {
        *self == Self::zero_mul()
    }
}

impl UnitMul for u32 {
    fn unit_mul() -> u32 {
        1
    }
    fn zero_mul() -> u32 {
        0
    }
}
impl UnitMul for usize {
    fn unit_mul() -> usize {
        1
    }
    fn zero_mul() -> usize {
        0
    }
}
impl UnitMul for f64 {
    fn unit_mul() -> f64 {
        1.0
    }
    fn zero_mul() -> f64 {
        0.0
    }
}
impl UnitMul for Prob {
    fn unit_mul() -> Prob {
        Prob::from_log_prob(0.0)
    }
    fn zero_mul() -> Prob {
        Prob::from_log_prob(f64::NEG_INFINITY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;

    #[test]
    fn unit_u32() {
        assert_eq!(u32::unit_add(), 0);
        assert_eq!(u32::unit_add() + 100, 100);
        assert_eq!(u32::unit_mul(), 1);
        assert_eq!(u32::unit_mul() * 100, 100);
        assert_eq!(u32::zero_mul(), 0);
        assert_eq!(u32::zero_mul() * 100, 0);
        assert!(!0u32.is_unit_mul());
        assert!(1u32.is_unit_mul());
        assert!(0u32.is_unit_add());
        assert!(0u32.is_zero_mul());
    }

    #[test]
    fn unit_f64() {
        assert_eq!(f64::unit_add(), 0.0);
        assert_eq!(f64::unit_add() + 100.0, 100.0);
        assert_eq!(f64::unit_mul(), 1.0);
        assert_eq!(f64::unit_mul() * 100.0, 100.0);
        assert_eq!(f64::zero_mul(), 0.0);
        assert_eq!(f64::zero_mul() * 100.0, 0.0);
        assert!(!0f64.is_unit_mul());
        assert!(!0.001.is_unit_add());
        assert!(1f64.is_unit_mul());
        assert!(0f64.is_unit_add());
        assert!(0f64.is_zero_mul());
    }

    #[test]
    fn unit_prob() {
        assert_eq!(Prob::unit_add(), p(0.0));
        assert_eq!(Prob::unit_add() + p(0.5), p(0.5));
        assert_eq!(Prob::unit_mul(), p(1.0));
        assert_eq!(Prob::unit_mul() * p(0.5), p(0.5));
        assert_eq!(Prob::zero_mul(), p(0.0));
        assert_eq!(Prob::zero_mul() * p(0.5), p(0.0));
        assert_eq!(p(0.0), p(0.0));
        assert_eq!(p(1.0), p(1.0));
        assert_eq!(p(0.5), p(0.5));
        assert!(p(1.0).is_unit_mul());
        assert!(!p(0.99999).is_unit_mul());
        assert!(p(0.0).is_unit_add());
        assert!(!p(0.0000001).is_unit_add());
        assert!(p(0.0).is_zero_mul());
        assert!(!p(0.0000001).is_zero_mul());
    }
}

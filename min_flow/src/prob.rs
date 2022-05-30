///
/// probability calculation
/// implements logaddexp
///
use approx::AbsDiffEq;

///
/// Wrapper of f64 that represents probability `0 <= p <= 1`
///
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Prob(f64);

///
/// short-hand of `Prob::from_prob`
///
pub fn p(p: f64) -> Prob {
    Prob::from_prob(p)
}

///
/// short-hand of `Prob::from_log_prob`
///
pub fn lp(lp: f64) -> Prob {
    Prob::from_log_prob(lp)
}

impl Prob {
    pub fn from_prob(value: f64) -> Prob {
        Prob(value.ln())
    }
    pub fn from_log_prob(log_value: f64) -> Prob {
        Prob(log_value)
    }
    ///
    /// Get the probability (in `[0, 1]`)
    pub fn to_value(self) -> f64 {
        self.0.exp()
    }
    ///
    /// Get the log probability
    pub fn to_log_value(self) -> f64 {
        self.0
    }
    ///
    /// Is `p == 0` or not?
    pub fn is_zero(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }
    ///
    /// prob=0.0
    ///
    pub fn zero() -> Prob {
        Prob(f64::NEG_INFINITY)
    }
    ///
    /// prob=1.0
    ///
    pub fn one() -> Prob {
        Prob(0.0)
    }
}

/// p=0 (Prob(-inf)) as a default value
impl Default for Prob {
    fn default() -> Self {
        Prob(f64::NEG_INFINITY)
    }
}

// display
impl std::fmt::Display for Prob {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:.2}(=log({:.12}))", self.0, self.to_value())
    }
}

// operators
impl std::ops::Add for Prob {
    type Output = Self;
    /// Given x>y
    /// log(exp(x) + exp(y))
    /// = log(exp(x) (1 + exp(y-x)))
    /// = log(exp(x)) + log(1 + exp(y-x))
    /// = x + log(1 + exp(y-x))
    fn add(self, other: Self) -> Self {
        let x = self.0;
        let y = other.0;
        if x == y {
            Prob(x + 2f64.ln())
        } else if x > y {
            Prob(x + ((y - x).exp() + 1.0).ln())
        } else {
            Prob(y + ((x - y).exp() + 1.0).ln())
        }
    }
}
impl std::ops::Mul for Prob {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Prob(self.0 + other.0)
    }
}
impl std::ops::Div for Prob {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Prob(self.0 - other.0)
    }
}
// assign
impl std::ops::AddAssign for Prob {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}
impl std::ops::MulAssign for Prob {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}
// sum/prod
impl std::iter::Sum for Prob {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(0.0), |a, b| a + b)
    }
}
impl<'a> std::iter::Sum<&'a Self> for Prob {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(0.0), |a, b| a + *b)
    }
}
impl std::iter::Product for Prob {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(1.0), |a, b| a * b)
    }
}
impl<'a> std::iter::Product<&'a Self> for Prob {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(1.0), |a, b| a * *b)
    }
}

/// for approx `assert_abs_diff_eq`
impl AbsDiffEq for Prob {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f64::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl Eq for Prob {}
impl Ord for Prob {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_id() {
//         let x = Prob::from_prob(0.3);
//         let e = Prob::from_prob(0.0);
//         assert_relative_eq!((x + e).0, x.0);
//         assert_relative_eq!((x * e).0, e.0);
//     }
//     #[test]
//     fn test_id_id() {
//         let x = Prob::from_prob(1.0);
//         let e = Prob::from_prob(0.0);
//         assert_relative_eq!((x + e).0, x.0);
//         assert_relative_eq!((x * e).0, e.0);
//     }
//     #[test]
//     fn test_sum() {
//         let xs = vec![
//             Prob::from_prob(0.1),
//             Prob::from_prob(0.1),
//             Prob::from_prob(0.1),
//         ];
//         let x: Prob = xs.iter().sum();
//         let y = Prob::from_prob(0.3);
//         assert_relative_eq!(x.to_value(), y.to_value());
//     }
//     #[test]
//     fn test_prod() {
//         let xs = vec![
//             Prob::from_prob(0.1),
//             Prob::from_prob(0.1),
//             Prob::from_prob(0.1),
//         ];
//         let x: Prob = xs.iter().product();
//         let y = Prob::from_prob(0.001);
//         assert_relative_eq!(x.to_value(), y.to_value());
//     }
//     #[test]
//     fn prob_add_mul() {
//         assert_eq!(p(0.0) + p(1.0), p(1.0));
//         assert_eq!(p(0.0) * p(1.0), p(0.0));
//         assert_abs_diff_eq!((p(0.3) + p(0.3)).0, p(0.6).0);
//         assert_abs_diff_eq!((p(0.3) * p(0.3)).0, p(0.09).0);
//         assert_abs_diff_eq!((p(0.5) + p(0.00001)).0, p(0.50001).0);
//         assert_abs_diff_eq!((p(0.5) * p(0.00001)).0, p(0.000005).0);
//     }
//     #[test]
//     fn prob_sum_prod() {
//         // sum/prod of zero element vec
//         let xs: Vec<Prob> = vec![];
//         let sum: Prob = xs.iter().sum();
//         let product: Prob = xs.iter().product();
//         assert_eq!(sum, p(0.0));
//         assert_eq!(product, p(1.0));

//         // sum/prod of vec of p=0
//         let xs: Vec<Prob> = vec![p(0.0), p(0.0)];
//         let sum: Prob = xs.iter().sum();
//         let product: Prob = xs.iter().product();
//         assert_eq!(sum, p(0.0));
//         assert_eq!(product, p(0.0));
//     }
//     #[test]
//     fn test_reflect() {
//         let x = Prob::from_prob(1.0);
//         let e = Prob::from_prob(0.0);
//         assert_relative_eq!((x + e).0, x.0);
//         assert_relative_eq!((x * e).0, e.0);
//         assert_relative_eq!((e * e).0, e.0);
//     }
//     #[test]
//     fn test_zero() {
//         let zero = Prob::from_prob(0.0);
//         println!("{:?}", zero);
//         assert!(zero.is_zero());
//         let nonzero = Prob::from_prob(0.00001);
//         assert!(!nonzero.is_zero());
//     }
//     #[test]
//     fn test_prob_assign() {
//         let mut x = p(0.4);
//         let y = p(0.2);
//         x += y;
//         assert_abs_diff_eq!(x, p(0.6));
//         let z = p(0.5);
//         x *= z;
//         assert_abs_diff_eq!(x, p(0.3));
//         let o = p(1.0);
//         x *= o;
//         assert_abs_diff_eq!(x, p(0.3));
//         let z = p(0.0);
//         x += z;
//         assert_abs_diff_eq!(x, p(0.3));
//         x *= z;
//         assert!(x.is_zero());
//     }
//     #[test]
//     fn prob_sort() {
//         // Sort by Ord and Eq
//         let mut ps = vec![p(0.9), p(0.2), p(0.5), p(0.1), p(1.0), p(0.0)];
//         ps.sort();
//         println!("{:?}", ps);
//         assert_eq!(ps[0], p(0.0));
//         assert_eq!(ps[1], p(0.1));
//         assert_eq!(ps[2], p(0.2));
//         assert_eq!(ps[3], p(0.5));
//         assert_eq!(ps[4], p(0.9));
//         assert_eq!(ps[5], p(1.0));
//     }
//     #[test]
//     fn prob_max_min() {
//         let mut ps = vec![p(0.9), p(0.2), p(0.5), p(0.1), p(1.0), p(0.0)];
//         let max = ps.iter().max().unwrap();
//         assert_eq!(*max, p(1.0));
//         let min = ps.iter().min().unwrap();
//         assert_eq!(*min, p(0.0));

//         assert!(p(0.1) > p(0.09999));
//         assert!(p(0.1) < p(0.100001));
//         assert!(p(0.0) < p(0.01));
//         assert!(p(1.0) > p(0.01));
//     }
//     #[test]
//     fn prob_assert_eq() {
//         assert!(abs_diff_eq!(p(0.1), p(0.1)));
//         assert!(!abs_diff_eq!(p(0.1), p(0.2)));
//         assert!(!abs_diff_eq!(p(0.1), p(0.11)));
//         assert!(abs_diff_eq!(p(0.1), p(0.11), epsilon = 0.1));
//         assert!(abs_diff_eq!(p(1.0), p(1.0)));
//         // TODO
//         assert!(!abs_diff_eq!(p(0.0), p(0.0)));
//     }
//     #[test]
//     fn prob_zero_one() {
//         assert_eq!(Prob::one(), Prob::from_prob(1.0));
//         assert_eq!(Prob::zero(), Prob::from_prob(0.0));
//         assert!(Prob::zero().is_zero());
//     }
// }

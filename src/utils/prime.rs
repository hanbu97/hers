use num_bigint_dig::{prime::probably_prime, BigUint};

/// A trait for checking primality of integers
pub trait PrimeChecking {
    /// Returns whether the number is prime
    fn is_prime(&self) -> bool;
}

impl PrimeChecking for u64 {
    fn is_prime(&self) -> bool {
        probably_prime(&BigUint::from(*self), 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime() {
        assert!(2u64.is_prime());
        assert!(3u64.is_prime());
        assert!(5u64.is_prime());
        assert!(7u64.is_prime());
        assert!(17u64.is_prime());
        assert!(101u64.is_prime());
        assert!(4611686018326724609u64.is_prime());

        assert!(!0u64.is_prime());
        assert!(!1u64.is_prime());
        assert!(!4u64.is_prime());
        assert!(!6u64.is_prime());
        assert!(!100u64.is_prime());
        assert!(!4611686018326724607u64.is_prime());
    }
}

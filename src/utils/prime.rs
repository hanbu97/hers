use super::num_bigint_ext::probably_prime;
use num_bigint::BigUint;

/// A trait for checking primality of integers
pub trait PrimeChecking {
    /// Returns whether the number is prime
    fn is_prime(&self) -> bool;
}

impl PrimeChecking for rug::Integer {
    fn is_prime(&self) -> bool {
        self.is_probably_prime(25) != rug::integer::IsPrime::No
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_is_prime() {
        assert!(BigUint::from(0xffffffffffffffc5u64).is_prime());
        assert!(BigUint::from(18446744073709551629u128).is_prime());
        assert!(!BigUint::from(0xffffffffffffffffu64).is_prime());
    }
}

impl PrimeChecking for u32 {
    fn is_prime(&self) -> bool {
        probably_prime(&BigUint::from(*self), 0)
    }
}

impl PrimeChecking for u64 {
    fn is_prime(&self) -> bool {
        probably_prime(&BigUint::from(*self), 0)
    }
}

impl PrimeChecking for u128 {
    fn is_prime(&self) -> bool {
        probably_prime(&BigUint::from(*self), 0)
    }
}

impl PrimeChecking for BigUint {
    fn is_prime(&self) -> bool {
        probably_prime(self, 0)
    }
}

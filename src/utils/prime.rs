use super::num_bigint_ext::probably_prime;
use num_bigint::BigUint;

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

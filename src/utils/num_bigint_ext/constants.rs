use super::*;

/// u64_digit
pub const BITS: u64 = 64;

lazy_static::lazy_static! {
    pub(crate) static ref BIG_1: BigUint = BigUint::one();
    pub(crate) static ref BIG_2: BigUint = BigUint::from_u64(2).unwrap();
    pub(crate) static ref BIG_3: BigUint = BigUint::from_u64(3).unwrap();
    pub(crate) static ref BIG_64: BigUint = BigUint::from_u64(64).unwrap();
}

pub(crate) const PRIMES_A: u64 = 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 37;
pub(crate) const PRIMES_B: u64 = 29 * 31 * 41 * 43 * 47 * 53;

/// Records the primes < 64.
pub(crate) const PRIME_BIT_MASK: u64 = 1 << 2
    | 1 << 3
    | 1 << 5
    | 1 << 7
    | 1 << 11
    | 1 << 13
    | 1 << 17
    | 1 << 19
    | 1 << 23
    | 1 << 29
    | 1 << 31
    | 1 << 37
    | 1 << 41
    | 1 << 43
    | 1 << 47
    | 1 << 53
    | 1 << 59
    | 1 << 61;

/// u64_digit
pub const BITS: u64 = 64;
use rug::Integer;

lazy_static::lazy_static! {
    pub(crate) static ref BIG_1: Integer = Integer::from(1u64);
    pub(crate) static ref BIG_2: Integer = Integer::from(2u64);
    pub(crate) static ref BIG_3: Integer = Integer::from(3u64);
    pub(crate) static ref BIG_64: Integer = Integer::from(64u64);
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

use super::*;

use num_bigint::BigUint;
use num_traits::ToPrimitive;

/// Computes the constant for the BRed algorithm.
/// Returns [((2^128)/q)/(2^64), (2^128)/q mod 2^64].
pub fn compute_barrett_constants(q: u64) -> [u64; 2] {
    let barrett = ((BigUint::from(1u64) << 128usize) / BigUint::from(q))
        .to_u128()
        .unwrap();

    let mlo = barrett as u64;
    let mhi = (barrett >> 64) as u64;

    [mhi, mlo]
}

pub fn b_red_add(a: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (mhi, _) = mul_hi_lo(a, bred_constant[0]);
    let mut r = a.wrapping_sub(mhi.wrapping_mul(q));
    if r >= q {
        r = r.wrapping_sub(q);
    }
    r
}
/// Computes a mod q in constant time.
/// The result is between 0 and 2*q-1.
pub fn b_red_add_lazy(x: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (s0, _) = mul_hi_lo(x, bred_constant[0]);
    x.wrapping_sub(s0.wrapping_mul(q))
}

/// Computes x*y mod q.
pub fn b_red(x: u64, y: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (mhi, mlo) = mul_hi_lo(x, y);

    let mut r = mhi.wrapping_mul(bred_constant[0]);

    let (hhi, hlo) = mul_hi_lo(mlo, bred_constant[0]);
    r = r.wrapping_add(hhi);

    let (lhi, _) = mul_hi_lo(mlo, bred_constant[1]);

    let (s0, carry) = hlo.overflowing_add(lhi);
    r = r.wrapping_add(carry as u64);

    let (hhi, hlo) = mul_hi_lo(mhi, bred_constant[1]);
    r = r.wrapping_add(hhi);

    let (_, carry) = hlo.overflowing_add(s0);
    r = r.wrapping_add(carry as u64);

    let mut r = mlo.wrapping_sub(r.wrapping_mul(q));

    if r >= q {
        r = r.wrapping_sub(q);
    }

    r
}

/// Computes x*y mod q in constant time.
/// The result is between 0 and 2*q-1.
pub fn b_red_lazy(x: u64, y: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (mhi, mlo) = mul_hi_lo(x, y);

    let mut r = mhi.wrapping_mul(bred_constant[0]);

    let (hhi, hlo) = mul_hi_lo(mlo, bred_constant[0]);
    r = r.wrapping_add(hhi);

    let (lhi, _) = mul_hi_lo(mlo, bred_constant[1]);

    let (s0, carry) = hlo.overflowing_add(lhi);
    r = r.wrapping_add(carry as u64);

    let (hhi, hlo) = mul_hi_lo(mhi, bred_constant[1]);
    r = r.wrapping_add(hhi);

    let (_, carry) = hlo.overflowing_add(s0);
    r = r.wrapping_add(carry as u64);

    mlo.wrapping_sub(r.wrapping_mul(q))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_b_red_add() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0),
            (1, 1),
            (2, 2),
            (18446744069414584320, 18446744069414584320),
            (18446744069414584321, 0),
            (18446744069414584322, 1),
            (4294967296, 4294967296),
            (9223372036854775808, 9223372036854775808),
            (18446744073709551615, 4294967294),
        ];

        for &(input, expected) in &test_cases {
            let result = b_red_add(input, q, bred_constant);
            assert_eq!(
                result, expected,
                "b_red_add({}) = {}; expected {}",
                input, result, expected
            );
        }
    }

    #[test]
    fn test_b_red_add_lazy() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0),
            (1, 1),
            (2, 2),
            (18446744069414584320, 18446744069414584320),
            (18446744069414584321, 18446744069414584321),
            (18446744069414584322, 18446744069414584322),
            (4294967296, 4294967296),
            (9223372036854775808, 9223372036854775808),
            (18446744073709551615, 18446744073709551615),
        ];

        for &(input, expected) in &test_cases {
            let result = b_red_add_lazy(input, q, bred_constant);
            assert_eq!(
                result, expected,
                "b_red_add_lazy({}) = {}; expected {}",
                input, result, expected
            );
        }
    }

    #[test]
    fn test_b_red() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0, 0),
            (1, 1, 1),
            (2, 2, 4),
            (18446744069414584320, 18446744069414584320, 1),
            (18446744069414584321, 18446744069414584321, 0),
            (18446744069414584322, 18446744069414584322, 1),
            (4294967296, 4294967296, 4294967295),
            (
                9223372036854775808,
                9223372036854775808,
                18446744068340842497,
            ),
            (
                18446744073709551615,
                18446744073709551615,
                18446744052234715141,
            ),
        ];

        for &(x, y, expected) in &test_cases {
            let result = b_red(x, y, q, bred_constant);
            assert_eq!(
                result, expected,
                "b_red({}, {}) = {}; expected {}",
                x, y, result, expected
            );
        }
    }

    #[test]
    fn test_b_red_lazy() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0, 0),
            (1, 1, 1),
            (2, 2, 4),
            (
                18446744069414584320,
                18446744069414584320,
                18446744069414584322,
            ),
            (
                18446744069414584321,
                18446744069414584321,
                18446744069414584321,
            ),
            (
                18446744069414584322,
                18446744069414584322,
                18446744069414584322,
            ),
            (4294967296, 4294967296, 4294967295),
            (
                9223372036854775808,
                9223372036854775808,
                18446744068340842497,
            ),
            (
                18446744073709551615,
                18446744073709551615,
                18446744052234715141,
            ),
        ];

        for &(x, y, expected) in &test_cases {
            let result = b_red_lazy(x, y, q, bred_constant);
            assert_eq!(
                result, expected,
                "b_red_lazy({}, {}) = {}; expected {}",
                x, y, result, expected
            );
        }
    }

    #[test]
    fn test_compute_barrett_constants() {
        let test_cases = [
            (3, 6148914691236517205, 6148914691236517205),
            (17, 1085102592571150095, 1085102592571150095),
            (257, 71777214294589695, 71777214294589695),
            (65537, 281470681808895, 281470681808895),
            (4294967295, 4294967297, 4294967297),
            (4294967296, 4294967296, 0),
            (9223372036854775807, 2, 4),
            (18446744073709551615, 1, 1),
            (0xFFFFFFFF00000001, 1, 4294967295),
        ];

        for &(q, expected_hi, expected_lo) in &test_cases {
            let constants = compute_barrett_constants(q);
            // println!("q: {}", q);
            // println!("((2^128)/q)/(2^64): {}", constants[0]);
            // println!("(2^128)/q mod 2^64: {}\n", constants[1]);

            assert_eq!(
                constants[0], expected_hi,
                "Mismatch in high 64 bits for q = {}",
                q
            );
            assert_eq!(
                constants[1], expected_lo,
                "Mismatch in low 64 bits for q = {}",
                q
            );
        }
    }
}

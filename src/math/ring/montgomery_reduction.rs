use super::helpers::mul_hi_lo;
use num_traits::WrappingMul;
use std::num::Wrapping;

/// Computes the constant m_red_constant = (q^-1) mod 2^64 required for MRed.
pub fn compute_montgomery_constant(q: u64) -> u64 {
    let mut mredconstant = 1u64;
    let mut q_temp = q;
    for _ in 0..63 {
        mredconstant = mredconstant.wrapping_mul(q_temp);
        q_temp = q_temp.wrapping_mul(q_temp);
    }
    mredconstant
}

/// Switches `a` to the Montgomery domain by computing a*2^64 mod q.
#[inline(always)]
pub fn m_form(a: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (hi, _lo) = mul_hi_lo(a, bred_constant[1]);
    let r = Wrapping(a)
        .wrapping_mul(&Wrapping(bred_constant[0]))
        .0
        .wrapping_add(hi)
        .wrapping_neg()
        .wrapping_mul(q);

    if r >= q {
        r.wrapping_sub(q)
    } else {
        r
    }
}

/// Switches `a` to the Montgomery domain by computing a*2^64 mod q in constant time.
/// The result is between 0 and 2*q-1.
#[inline(always)]
pub fn m_form_lazy(a: u64, q: u64, bred_constant: [u64; 2]) -> u64 {
    let (hi, _lo) = mul_hi_lo(a, bred_constant[1]);
    Wrapping(a)
        .wrapping_mul(&Wrapping(bred_constant[0]))
        .0
        .wrapping_add(hi)
        .wrapping_neg()
        .wrapping_mul(q)
}

/// Switches a from the Montgomery domain back to the standard domain by computing a*(1/2^64) mod q.
pub fn i_m_form(a: u64, q: u64, mred_constant: u64) -> u64 {
    let (r, _) = mul_hi_lo(a.wrapping_mul(mred_constant), q);
    let mut result = q.wrapping_sub(r);
    if result >= q {
        result = result.wrapping_sub(q);
    }
    result
}

/// Switches a from the Montgomery domain back to the standard domain by computing a*(1/2^64) mod q in constant time.
/// The result is between 0 and 2*q-1.
pub fn i_m_form_lazy(a: u64, q: u64, mred_constant: u64) -> u64 {
    let (r, _) = mul_hi_lo(a.wrapping_mul(mred_constant), q);
    q.wrapping_sub(r)
}

#[test]
fn test_compute_montgomery_constant() {
    let test_cases = [
        (3, 12297829382473034411),
        (17, 17361641481138401521),
        (257, 18374966859414961921),
        (65537, 18446462603027742721),
        (4294967295, 18446744069414584319),           // 2^32 - 1
        (4294967296, 0),                              // 2^32
        (9223372036854775807, 9223372036854775807),   // 2^63 - 1
        (18446744073709551615, 18446744073709551615), // 2^64 - 1
    ];

    for &(q, expected_m_red_constant) in &test_cases {
        let m_red_constant = compute_montgomery_constant(q);

        // println!("q: {}", q);
        // println!("m_red_constant: {}", m_red_constant);
        // println!("expected: {}\n", expected_m_red_constant);

        assert_eq!(
            m_red_constant, expected_m_red_constant,
            "Mismatch for q = {}: expected {}, but got {}",
            q, expected_m_red_constant, m_red_constant
        );
    }
}

/// Computes x * y * (1/2^64) mod q.
pub fn m_red(x: u64, y: u64, q: u64, mred_constant: u64) -> u64 {
    let (mhi, mlo) = mul_hi_lo(x, y);
    let (hhi, _) = mul_hi_lo(mlo.wrapping_mul(mred_constant), q);
    let mut r = mhi.wrapping_sub(hhi).wrapping_add(q);
    if r >= q {
        r = r.wrapping_sub(q);
    }
    r
}

/// Computes x * y * (1/2^64) mod q in constant time.
/// The result is between 0 and 2*q-1.
pub fn m_red_lazy(x: u64, y: u64, q: u64, mred_constant: u64) -> u64 {
    let (ahi, alo) = mul_hi_lo(x, y);
    let (h, _) = mul_hi_lo(alo.wrapping_mul(mred_constant), q);
    ahi.wrapping_sub(h).wrapping_add(q)
}

#[cfg(test)]
mod tests {
    use crate::math::ring::barrett_reduction::compute_barrett_constants;

    use super::*;

    #[test]
    fn test_m_red_lazy() {
        let q = 0xFFFFFFFF00000001;
        let mred_constant = compute_montgomery_constant(q);

        let test_cases = [
            (0, 0, 18446744069414584321),
            (0, 1, 18446744069414584321),
            (1, 1, 18446744065119617025),
            (2, 2, 18446744052234715137),
            (3, 3, 18446744030759878657),
            (4, 4, 18446744000695107585),
            (5, 5, 18446743962040401921),
            (
                18446744069414584319,
                18446744069414584319,
                18446744047939747842,
            ),
            (
                18446744069414584320,
                18446744069414584320,
                18446744060824649730,
            ),
            (
                18446744069414584321,
                18446744069414584321,
                18446744069414584321,
            ),
            (
                18446744069414584322,
                18446744069414584322,
                18446744060824649730,
            ),
            (
                18446744069414584323,
                18446744069414584323,
                18446744047939747842,
            ),
            (4294967295, 4294967295, 4294967295),
            (4294967296, 4294967296, 18446744069414584322),
            (4294967297, 4294967297, 18446744056529682436),
            (
                9223372036854775807,
                9223372036854775807,
                4611686009837453312,
            ),
            (
                9223372036854775808,
                9223372036854775808,
                4611686014132420609,
            ),
            (
                9223372036854775809,
                9223372036854775809,
                4611686009837453314,
            ),
            (
                18446744073709551615,
                18446744073709551615,
                18446744065119617023,
            ),
        ];

        for &(x, y, expected) in &test_cases {
            let result = m_red_lazy(x, y, q, mred_constant);
            assert_eq!(
                result, expected,
                "m_red_lazy({}, {}) = {}; expected {}",
                x, y, result, expected
            );
        }
    }

    #[test]
    fn test_m_red() {
        let q = 0xFFFFFFFF00000001;
        let mred_constant = compute_montgomery_constant(q);

        let test_cases = [
            (0, 0, 0),
            (0, 1, 0),
            (1, 1, 18446744065119617025),
            (2, 2, 18446744052234715137),
            (3, 3, 18446744030759878657),
            (4, 4, 18446744000695107585),
            (5, 5, 18446743962040401921),
            (
                18446744069414584319,
                18446744069414584319,
                18446744047939747842,
            ),
            (
                18446744069414584320,
                18446744069414584320,
                18446744060824649730,
            ),
            (18446744069414584321, 18446744069414584321, 0),
            (
                18446744069414584322,
                18446744069414584322,
                18446744060824649730,
            ),
            (
                18446744069414584323,
                18446744069414584323,
                18446744047939747842,
            ),
            (4294967295, 4294967295, 4294967295),
            (4294967296, 4294967296, 1),
            (4294967297, 4294967297, 18446744056529682436),
            (
                9223372036854775807,
                9223372036854775807,
                4611686009837453312,
            ),
            (
                9223372036854775808,
                9223372036854775808,
                4611686014132420609,
            ),
            (
                9223372036854775809,
                9223372036854775809,
                4611686009837453314,
            ),
            (
                18446744073709551615,
                18446744073709551615,
                18446744065119617023,
            ),
        ];

        for &(x, y, expected) in &test_cases {
            let result = m_red(x, y, q, mred_constant);
            assert_eq!(
                result, expected,
                "m_red({}, {}) = {}; expected {}",
                x, y, result, expected
            );
        }
    }

    #[test]
    fn test_m_form() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0),
            (1, 4294967295),
            (2, 8589934590),
            (3, 12884901885),
            (4, 17179869180),
            (5, 21474836475),
            (18446744069414584319, 18446744060824649731),
            (18446744069414584320, 18446744065119617026),
            (18446744069414584321, 0),
            (18446744069414584322, 0),
            (18446744069414584323, 4294967295),
            (4294967295, 18446744065119617025),
            (4294967296, 18446744069414584320),
            (4294967297, 4294967294),
            (9223372036854775807, 18446744062972133378),
            (9223372036854775808, 18446744067267100673),
            (9223372036854775809, 2147483647),
            (18446744073709551615, 18446744056529682435),
        ];

        println!("Testing m_form:");
        for &(input, expected) in &test_cases {
            let result = m_form(input, q, bred_constant);
            println!("m_form({}) = {}", input, result);
            assert_eq!(
                result, expected,
                "m_form({}) = {}; expected {}",
                input, result, expected
            );
        }
    }

    #[test]
    fn test_m_form_lazy() {
        let q = 0xFFFFFFFF00000001;
        let bred_constant = compute_barrett_constants(q);

        let test_cases = [
            (0, 0),
            (1, 4294967295),
            (2, 8589934590),
            (3, 12884901885),
            (4, 17179869180),
            (5, 21474836475),
            (18446744069414584319, 18446744060824649731),
            (18446744069414584320, 18446744065119617026),
            (18446744069414584321, 18446744069414584321),
            (18446744069414584322, 0),
            (18446744069414584323, 4294967295),
            (4294967295, 18446744065119617025),
            (4294967296, 18446744069414584320),
            (4294967297, 18446744073709551615),
            (9223372036854775807, 18446744062972133378),
            (9223372036854775808, 18446744067267100673),
            (9223372036854775809, 18446744071562067968),
            (18446744073709551615, 18446744056529682435),
        ];

        println!("\nTesting m_form_lazy:");
        for &(input, expected) in &test_cases {
            let result = m_form_lazy(input, q, bred_constant);
            println!("m_form_lazy({}) = {}", input, result);
            assert_eq!(
                result, expected,
                "m_form_lazy({}) = {}; expected {}",
                input, result, expected
            );
        }
    }

    #[test]
    fn test_i_m_form() {
        let q = 0xFFFFFFFF00000001;
        let mred_constant = compute_montgomery_constant(q);

        let test_cases = [
            (0, 0),
            (1, 18446744065119617025),
            (2, 18446744060824649729),
            (3, 18446744056529682433),
            (4, 18446744052234715137),
            (5, 18446744047939747841),
            (18446744069414584319, 8589934592),
            (18446744069414584320, 4294967296),
            (18446744069414584321, 0),
            (18446744069414584322, 18446744065119617025),
            (18446744069414584323, 18446744060824649729),
            (4294967295, 1),
            (4294967296, 18446744065119617026),
            (4294967297, 18446744060824649730),
            (9223372036854775807, 9223372039002259457),
            (9223372036854775808, 9223372034707292161),
            (9223372036854775809, 9223372030412324865),
            (18446744073709551615, 4294967297),
        ];

        println!("Testing i_m_form:");
        for &(input, expected) in &test_cases {
            let result = i_m_form(input, q, mred_constant);
            println!("i_m_form({}) = {}", input, result);
            assert_eq!(
                result, expected,
                "i_m_form({}) = {}; expected {}",
                input, result, expected
            );
        }
    }

    #[test]
    fn test_i_m_form_lazy() {
        let q = 0xFFFFFFFF00000001;
        let mred_constant = compute_montgomery_constant(q);

        let test_cases = [
            (0, 18446744069414584321),
            (1, 18446744065119617025),
            (2, 18446744060824649729),
            (3, 18446744056529682433),
            (4, 18446744052234715137),
            (5, 18446744047939747841),
            (18446744069414584319, 8589934592),
            (18446744069414584320, 4294967296),
            (18446744069414584321, 18446744069414584321),
            (18446744069414584322, 18446744065119617025),
            (18446744069414584323, 18446744060824649729),
            (4294967295, 1),
            (4294967296, 18446744065119617026),
            (4294967297, 18446744060824649730),
            (9223372036854775807, 9223372039002259457),
            (9223372036854775808, 9223372034707292161),
            (9223372036854775809, 9223372030412324865),
            (18446744073709551615, 4294967297),
        ];

        println!("\nTesting i_m_form_lazy:");
        for &(input, expected) in &test_cases {
            let result = i_m_form_lazy(input, q, mred_constant);
            println!("i_m_form_lazy({}) = {}", input, result);
            assert_eq!(
                result, expected,
                "i_m_form_lazy({}) = {}; expected {}",
                input, result, expected
            );
        }
    }
}

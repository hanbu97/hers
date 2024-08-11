use super::{
    reduction::{
        barrett::{b_red, b_red_add, b_red_add_lazy, b_red_lazy, compute_barrett_constants},
        conditional::c_red,
    },
    subring::SubRing,
};
use crate::math::ring::reduction::montgomery::*;

pub mod ntt_operations;
pub mod vec_operations;

/// Compute rescale constants for subring moduli
///
/// This function computes rescale constants for a given set of subrings.
/// These constants are used in the context of Residue Number System (RNS)
/// arithmetic, particularly for scaling between different moduli.
///
/// # Arguments
/// * `sub_rings` - A slice of SubRing structures representing different moduli
///
/// # Returns
/// A vector of vectors containing the computed rescale constants
pub fn compute_rescale_constants(sub_rings: &[SubRing]) -> Vec<Vec<u64>> {
    let mut rescale_constants = Vec::with_capacity(sub_rings.len() - 1);

    for j in (1..sub_rings.len()).rev() {
        let qj = sub_rings[j].modulus;
        let mut constants = Vec::with_capacity(j);

        for i in 0..j {
            let qi = sub_rings[i].modulus;
            let constant = m_form(
                qi.wrapping_sub(mod_exp(qj, qi.wrapping_sub(2), qi)),
                qi,
                sub_rings[i].b_red_constant,
            );
            constants.push(constant);
        }

        rescale_constants.push(constants);
    }

    rescale_constants.reverse();
    rescale_constants
}

/// Evaluates a polynomial at a given point x modulo p
///
/// # Arguments
/// * `x` - The point at which to evaluate the polynomial
/// * `poly` - The coefficients of the polynomial
/// * `p` - The modulus
///
/// # Returns
/// The result of the polynomial evaluation
pub fn eval_poly_mod_p(x: u64, poly: &[u64], p: u64) -> u64 {
    let brc = compute_barrett_constants(p);
    let mut y = poly[poly.len() - 1];
    for &coeff in poly.iter().rev().skip(1) {
        y = b_red(y, x, p, brc);
        y = c_red(y.wrapping_add(coeff), p);
    }
    y
}

/// Returns the minimum of two integers
///
/// # Arguments
/// * `x` - First integer
/// * `y` - Second integer
///
/// # Returns
/// The smaller of the two input integers
pub fn min(x: i32, y: i32) -> i32 {
    if x > y {
        y
    } else {
        x
    }
}

/// Performs modular exponentiation: x^e mod p
///
/// # Arguments
/// * `x` - Base
/// * `e` - Exponent
/// * `p` - Modulus
///
/// # Returns
/// The result of x^e mod p
pub fn mod_exp(x: u64, e: u64, p: u64) -> u64 {
    let brc = compute_barrett_constants(p);
    let mut result = 1;
    let mut base = x;
    let mut exp = e;
    while exp > 0 {
        if exp & 1 == 1 {
            result = b_red(result, base, p, brc);
        }
        base = b_red(base, base, p, brc);
        exp >>= 1;
    }
    result
}

/// Performs modular exponentiation when the modulus is a power of two
///
/// # Arguments
/// * `x` - Base
/// * `e` - Exponent
/// * `p` - Modulus (must be a power of two)
///
/// # Returns
/// The result of x^e mod p
pub fn mod_exp_pow2(x: u64, e: u64, p: u64) -> u64 {
    let mut result: u64 = 1;
    let mut base: u64 = x;
    let mut exp = e;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result.wrapping_mul(base);
        }
        base = base.wrapping_mul(base);
        exp >>= 1;
    }
    result & (p - 1)
}

/// Performs modular exponentiation using Montgomery reduction
///
/// # Arguments
/// * `x` - Base in Montgomery form
/// * `e` - Exponent
/// * `q` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `b_red_constant` - Barrett reduction constants
///
/// # Returns
/// The result of x^e mod q in Montgomery form
pub fn modexp_montgomery(
    x: u64,
    e: i32,
    q: u64,
    m_red_constant: u64,
    b_red_constant: [u64; 2],
) -> u64 {
    let mut result = m_form(1, q, b_red_constant);
    let mut base = x;
    let mut exp = e;
    while exp > 0 {
        if exp & 1 == 1 {
            result = m_red(result, base, q, m_red_constant);
        }
        base = m_red(base, base, q, m_red_constant);
        exp >>= 1;
    }
    result
}

use super::vec_operations::{
    mul_scalar_montgomery_lazy_vec_self, mul_scalar_montgomery_vec_self, reduce_vec_self,
};
use super::*;

/// Performs the NTT algorithm with final reduction
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `degree` - Ring size
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `b_red_constant` - Barrett reduction constants
/// * `roots` - NTT roots
pub fn ntt_standard(
    p1: &[u64],
    p2: &mut [u64],
    degree: u64,
    modulus: u64,
    m_red_constant: u64,
    b_red_constant: [u64; 2],
    roots: &[u64],
) {
    println!("ntt_standard roots: {:?}", roots);
    ntt_core_lazy(p1, p2, degree as usize, modulus, m_red_constant, roots);
    reduce_vec_self(p2, modulus, b_red_constant);
}

/// Performs the NTT algorithm without final reduction
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - NTT roots
pub fn ntt_standard_lazy(
    p1: &[u64],
    p2: &mut [u64],
    n: usize,
    modulus: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    println!("ntt_standard_lazy");
    ntt_core_lazy(p1, p2, n, modulus, m_red_constant, roots);
}

/// Performs the INTT algorithm with final reduction and multiplication by N^-1
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `n_inv` - N^-1 mod q in Montgomery form
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - INTT roots
pub fn intt_standard(
    p1: &[u64],
    p2: &mut [u64],
    n: usize,
    n_inv: u64,
    modulus: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    intt_core_lazy(p1, p2, n, modulus, m_red_constant, roots);
    mul_scalar_montgomery_vec_self(p2, n_inv, modulus, m_red_constant);
}

/// Performs the core NTT algorithm with lazy reduction
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - NTT roots
pub fn ntt_core_lazy(
    p1: &[u64],
    p2: &mut [u64],
    degree: usize,
    modulus: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    if p1.len() < degree || p2.len() < degree || roots.len() < degree {
        panic!(
            "cannot ntt_core_lazy: ensure that len(p1)={}, len(p2)={} and len(roots)={} >= N={}",
            p1.len(),
            p2.len(),
            roots.len(),
            degree
        );
    }

    ntt_lazy(p1, p2, degree, modulus, m_red_constant, roots);
}

/// Performs the INTT algorithm without final reduction and multiplication by N^-1
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `n_inv` - N^-1 mod q in Montgomery form
/// * `q` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - INTT roots
pub fn intt_standard_lazy(
    p1: &[u64],
    p2: &mut [u64],
    n: usize,
    n_inv: u64,
    q: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    intt_core_lazy(p1, p2, n, q, m_red_constant, roots);
    mul_scalar_montgomery_lazy_vec_self(p2, n_inv, q, m_red_constant);
}

/// Performs the core INTT algorithm with lazy reduction
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `q` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - INTT roots
pub fn intt_core_lazy(
    p1: &[u64],
    p2: &mut [u64],
    n: usize,
    q: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    if p1.len() < n || p2.len() < n || roots.len() < n {
        panic!(
            "cannot intt_core_lazy: ensure that len(p1)={}, len(p2)={} and len(roots)={} >= N={}",
            p1.len(),
            p2.len(),
            roots.len(),
            n
        );
    }

    intt_lazy(p1, p2, n, q, m_red_constant, roots);
}

/// Performs the butterfly operation in the NTT algorithm
///
/// Formula: X = U + V * Psi mod Q, Y = U - V * Psi mod Q
///
/// # Arguments
/// * `u`, `v` - Input values
/// * `psi` - Twiddle factor
/// * `two_q` - 2 * Q
/// * `four_q` - 4 * Q
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
pub fn butterfly(
    u: u64,
    v: u64,
    psi: u64,
    two_q: u64,
    four_q: u64,
    modulus: u64,
    m_red_constant: u64,
) -> (u64, u64) {
    let u = if u >= four_q { u - four_q } else { u };
    let v = m_red_lazy(v, psi, modulus, m_red_constant);
    (u + v, u + two_q - v)
}

/// Performs the inverse butterfly operation in the INTT algorithm
///
/// Formula: X = (U + V) mod Q, Y = (U - V) * Psi mod Q
///
/// # Arguments
/// * `u`, `v` - Input values
/// * `psi` - Twiddle factor
/// * `two_q` - 2 * Q
/// * `four_q` - 4 * Q
/// * `q` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
pub fn inv_butterfly(
    u: u64,
    v: u64,
    psi: u64,
    two_q: u64,
    four_q: u64,
    q: u64,
    m_red_constant: u64,
) -> (u64, u64) {
    let x = u + v;
    let x = if x >= two_q { x - two_q } else { x };
    let y = m_red_lazy(u + four_q - v, psi, q, m_red_constant);
    (x, y)
}

/// Performs the NTT algorithm with lazy reduction
///
/// This function implements the Cooley-Tukey FFT algorithm for NTT.
/// The calculation is done in-place with the following steps:
/// 1. Initial butterfly operations
/// 2. Subsequent butterfly operations in a divide-and-conquer manner
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `modulus` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - NTT roots
pub fn ntt_lazy(
    p1: &[u64],
    p2: &mut [u64],
    n: usize,
    modulus: u64,
    m_red_constant: u64,
    roots: &[u64],
) {
    let mut t = n >> 1;
    let four_q = 4 * modulus;
    let two_q = 2 * modulus;

    let f = roots[1];
    for (jx, jy) in (0..t).zip(t..n) {
        let (x, y) = butterfly(p1[jx], p1[jy], f, two_q, four_q, modulus, m_red_constant);
        p2[jx] = x;
        p2[jy] = y;
    }

    let mut m = 2;
    while m < n {
        t >>= 1;

        for i in 0..m {
            let j1 = (i * t) << 1;
            let j2 = j1 + t;
            let f = roots[m + i];

            for (jx, jy) in (j1..j2).zip(j1 + t..j2 + t) {
                let (x, y) = butterfly(p2[jx], p2[jy], f, two_q, four_q, modulus, m_red_constant);
                p2[jx] = x;
                p2[jy] = y;
            }
        }

        m <<= 1;
    }
}

/// Performs the INTT algorithm with lazy reduction
///
/// This function implements the Gentleman-Sande FFT algorithm for INTT.
/// The calculation is done in-place with the following steps:
/// 1. Butterfly operations in a divide-and-conquer manner
/// 2. Final butterfly operations
///
/// # Arguments
/// * `p1` - Input slice
/// * `p2` - Output slice
/// * `n` - Ring size
/// * `q` - Modulus
/// * `m_red_constant` - Montgomery reduction constant
/// * `roots` - INTT roots
pub fn intt_lazy(p1: &[u64], p2: &mut [u64], n: usize, q: u64, m_red_constant: u64, roots: &[u64]) {
    let mut t = 1;
    let mut h = n >> 1;
    let two_q = q << 1;
    let four_q = q << 2;

    for i in 0..h {
        let f = roots[h + i];
        for (jx, jy) in (i * 2 * t..i * 2 * t + t).zip(i * 2 * t + t..i * 2 * t + 2 * t) {
            let (x, y) = inv_butterfly(p1[jx], p1[jy], f, two_q, four_q, q, m_red_constant);
            p2[jx] = x;
            p2[jy] = y;
        }
    }

    t <<= 1;

    let mut m = n >> 1;
    while m > 1 {
        h = m >> 1;

        for i in 0..h {
            let f = roots[h + i];
            for (jx, jy) in (i * 2 * t..i * 2 * t + t).zip(i * 2 * t + t..i * 2 * t + 2 * t) {
                let (x, y) = inv_butterfly(p2[jx], p2[jy], f, two_q, four_q, q, m_red_constant);
                p2[jx] = x;
                p2[jy] = y;
            }
        }

        t <<= 1;
        m >>= 1;
    }
}

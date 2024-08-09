use crate::math::ring::reduction::conditional::c_red;

use super::reduction::{
    barrett::{b_red, b_red_add, b_red_add_lazy, b_red_lazy},
    montgomery::{im_form, m_form, m_form_lazy, m_red, m_red_lazy},
};

#[inline(always)]
pub fn add_vec(p1: &[u64], p2: &[u64], p3: &mut [u64], modulus: u64) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| c_red(x.wrapping_add(y), modulus))
        .zip(p3.iter_mut())
        .for_each(|(sum, out)| *out = sum);
}

#[inline(always)]
pub fn add_lazy_vec(p1: &[u64], p2: &[u64], p3: &mut [u64]) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| x.wrapping_add(y))
        .zip(p3.iter_mut())
        .for_each(|(sum, out)| *out = sum);
}

#[inline(always)]
pub fn sub_vec(p1: &[u64], p2: &[u64], p3: &mut [u64], modulus: u64) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| c_red(x.wrapping_add(modulus).wrapping_sub(y), modulus))
        .zip(p3.iter_mut())
        .for_each(|(diff, out)| *out = diff);
}

#[inline(always)]
pub fn sub_lazy_vec(p1: &[u64], p2: &[u64], p3: &mut [u64], modulus: u64) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| x.wrapping_add(modulus).wrapping_sub(y))
        .zip(p3.iter_mut())
        .for_each(|(diff, out)| *out = diff);
}

#[inline(always)]
pub fn neg_vec(p1: &[u64], p2: &mut [u64], modulus: u64) {
    p1.iter()
        .map(|&x| modulus.wrapping_sub(x))
        .zip(p2.iter_mut())
        .for_each(|(neg, out)| *out = neg);
}

#[inline(always)]
pub fn reduce_vec(p1: &[u64], p2: &mut [u64], modulus: u64, bred_constant: [u64; 2]) {
    p1.iter()
        .map(|&x| b_red_add(x, modulus, bred_constant))
        .zip(p2.iter_mut())
        .for_each(|(reduced, out)| *out = reduced);
}

#[inline(always)]
pub fn reduce_lazy_vec(p1: &[u64], p2: &mut [u64], modulus: u64, bred_constant: [u64; 2]) {
    p1.iter()
        .map(|&x| b_red_add_lazy(x, modulus, bred_constant))
        .zip(p2.iter_mut())
        .for_each(|(reduced, out)| *out = reduced);
}

#[inline(always)]
pub fn mul_coeffs_lazy_vec(p1: &[u64], p2: &[u64], p3: &mut [u64]) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| x.wrapping_mul(y))
        .zip(p3.iter_mut())
        .for_each(|(product, out)| *out = product);
}

#[inline(always)]
pub fn mul_coeffs_lazy_then_add_lazy_vec(p1: &[u64], p2: &[u64], p3: &mut [u64]) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| *out = out.wrapping_add(x.wrapping_mul(y)));
}

#[inline(always)]
pub fn mul_coeffs_barrett_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    bred_constant: [u64; 2],
) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| b_red(x, y, modulus, bred_constant))
        .zip(p3.iter_mut())
        .for_each(|(product, out)| *out = product);
}

#[inline(always)]
pub fn mul_coeffs_barrett_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    bred_constant: [u64; 2],
) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| b_red_lazy(x, y, modulus, bred_constant))
        .zip(p3.iter_mut())
        .for_each(|(product, out)| *out = product);
}

#[inline(always)]
pub fn mul_coeffs_then_add_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    bred_constant: [u64; 2],
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| {
            *out = c_red(
                out.wrapping_add(b_red(x, y, modulus, bred_constant)),
                modulus,
            )
        });
}

#[inline(always)]
pub fn mul_coeffs_barrett_then_add_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    bred_constant: [u64; 2],
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| *out = out.wrapping_add(b_red(x, y, modulus, bred_constant)));
}

#[inline(always)]
pub fn mul_coeffs_montgomery_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| m_red(x, y, modulus, mred_constant))
        .zip(p3.iter_mut())
        .for_each(|(product, out)| *out = product);
}

#[inline(always)]
pub fn mul_coeffs_montgomery_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .map(|(&x, &y)| m_red_lazy(x, y, modulus, mred_constant))
        .zip(p3.iter_mut())
        .for_each(|(product, out)| *out = product);
}

#[inline(always)]
pub fn mul_coeffs_montgomery_then_add_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| {
            *out = c_red(
                out.wrapping_add(m_red(x, y, modulus, mred_constant)),
                modulus,
            )
        });
}

#[inline(always)]
pub fn mul_coeffs_montgomery_then_add_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| *out = out.wrapping_add(m_red(x, y, modulus, mred_constant)));
}

#[inline(always)]
pub fn mul_coeffs_montgomery_lazy_then_add_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), out)| {
            *out = out.wrapping_add(m_red_lazy(x, y, modulus, mred_constant))
        });
}

#[inline(always)]
pub fn mul_coeffs_montgomery_then_sub_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = c_red(
                z.wrapping_add(modulus.wrapping_sub(m_red(x, y, modulus, mred_constant))),
                modulus,
            );
        });
}

#[inline(always)]
pub fn mul_coeffs_montgomery_then_sub_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = z.wrapping_add(modulus.wrapping_sub(m_red(x, y, modulus, mred_constant)));
        });
}

#[inline(always)]
pub fn mul_coeffs_montgomery_lazy_then_sub_lazy_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    let two_modulus = modulus << 1;
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = z.wrapping_add(two_modulus.wrapping_sub(m_red_lazy(x, y, modulus, mred_constant)));
        });
}

#[inline(always)]
pub fn mul_coeffs_montgomery_lazy_then_neg_vec(
    p1: &[u64],
    p2: &[u64],
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    let two_modulus = modulus << 1;
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = two_modulus.wrapping_sub(m_red_lazy(x, y, modulus, mred_constant));
        });
}

#[inline(always)]
pub fn add_lazy_then_mul_scalar_montgomery_vec(
    p1: &[u64],
    p2: &[u64],
    scalar_mont: u64,
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = m_red(x.wrapping_add(y), scalar_mont, modulus, mred_constant);
        });
}

#[inline(always)]
pub fn add_scalar_lazy_then_mul_scalar_montgomery_vec(
    p1: &[u64],
    scalar0: u64,
    scalar_mont1: u64,
    p2: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = m_red(
            x.wrapping_add(scalar0),
            scalar_mont1,
            modulus,
            mred_constant,
        );
    });
}

#[inline(always)]
pub fn add_scalar_vec(p1: &[u64], scalar: u64, p2: &mut [u64], modulus: u64) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = c_red(x.wrapping_add(scalar), modulus);
    });
}

#[inline(always)]
pub fn add_scalar_lazy_vec(p1: &[u64], scalar: u64, p2: &mut [u64]) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = x.wrapping_add(scalar);
    });
}

#[inline(always)]
pub fn add_scalar_lazy_then_neg_two_modulus_lazy_vec(
    p1: &[u64],
    scalar: u64,
    p2: &mut [u64],
    modulus: u64,
) {
    let two_modulus = modulus << 1;
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = scalar.wrapping_add(two_modulus).wrapping_sub(x);
    });
}

#[inline(always)]
pub fn sub_scalar_vec(p1: &[u64], scalar: u64, p2: &mut [u64], modulus: u64) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = c_red(x.wrapping_add(modulus).wrapping_sub(scalar), modulus);
    });
}

#[inline(always)]
pub fn mul_scalar_montgomery_vec(
    p1: &[u64],
    scalar_mont: u64,
    p2: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = m_red(x, scalar_mont, modulus, mred_constant);
    });
}

#[inline(always)]
pub fn mul_scalar_montgomery_lazy_vec(
    p1: &[u64],
    scalar_mont: u64,
    p2: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = m_red_lazy(x, scalar_mont, modulus, mred_constant);
    });
}

#[inline(always)]
pub fn mul_scalar_montgomery_then_add_vec(
    p1: &[u64],
    scalar_mont: u64,
    p2: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = c_red(
            z.wrapping_add(m_red(x, scalar_mont, modulus, mred_constant)),
            modulus,
        );
    });
}

#[inline(always)]
pub fn mul_scalar_montgomery_then_add_scalar_vec(
    p1: &[u64],
    scalar0: u64,
    scalar_mont1: u64,
    p2: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = c_red(
            m_red(x, scalar_mont1, modulus, mred_constant).wrapping_add(scalar0),
            modulus,
        );
    });
}

#[inline(always)]
pub fn sub_then_mul_scalar_montgomery_two_modulus_vec(
    p1: &[u64],
    p2: &[u64],
    scalar_mont: u64,
    p3: &mut [u64],
    modulus: u64,
    mred_constant: u64,
) {
    let two_modulus = modulus << 1;
    p1.iter()
        .zip(p2)
        .zip(p3.iter_mut())
        .for_each(|((&x, &y), z)| {
            *z = m_red(
                two_modulus.wrapping_sub(y).wrapping_add(x),
                scalar_mont,
                modulus,
                mred_constant,
            );
        });
}

#[inline(always)]
pub fn m_form_vec(p1: &[u64], p2: &mut [u64], modulus: u64, bred_constant: [u64; 2]) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = m_form(x, modulus, bred_constant);
    });
}

#[inline(always)]
pub fn m_form_lazy_vec(p1: &[u64], p2: &mut [u64], modulus: u64, bred_constant: [u64; 2]) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = m_form_lazy(x, modulus, bred_constant);
    });
}

#[inline(always)]
pub fn im_form_vec(p1: &[u64], p2: &mut [u64], modulus: u64, mred_constant: u64) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = im_form(x, modulus, mred_constant);
    });
}

#[inline(always)]
pub fn zero_vec(p1: &mut [u64]) {
    p1.iter_mut().for_each(|z| *z = 0);
}

#[inline(always)]
pub fn mask_vec(p1: &[u64], w: u32, mask: u64, p2: &mut [u64]) {
    p1.iter().zip(p2.iter_mut()).for_each(|(&x, z)| {
        *z = (x >> w) & mask;
    });
}

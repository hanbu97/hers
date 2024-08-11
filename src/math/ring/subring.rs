use super::{errors::SubRingError, operations::vec_operations::*};
use crate::math::{
    ntt::{
        params::{NTTParams, NTTTable},
        traits::NumberTheoreticTransform,
        NTTImplementations,
    },
    ring::{
        constants::MINIMUM_RING_DEGREE_FOR_LOOP_UNROLLED_OPS,
        reduction::{barrett::compute_barrett_constants, montgomery::compute_montgomery_constant},
    },
};

fn calculate_mask(modulus: u64) -> u64 {
    (1u64 << (64 - (modulus - 1).leading_zeros())) - 1
}

/// SubRing stores precomputation for fast modular reduction and NTT for a given modulus.
/// It encapsulates all the necessary information and methods to perform efficient
/// polynomial operations in a specific prime field.
pub struct SubRing {
    /// The Number Theoretic Transform implementation.
    /// This provides flexibility in choosing different NTT algorithms.
    pub ntt: NTTImplementations,

    /// Number of coefficients in the polynomial.
    /// This is typically a power of 2 for efficient NTT operations.
    pub n: usize,

    /// The prime modulus defining the finite field.
    pub modulus: u64,

    /// Unique prime factors of (modulus - 1).
    /// Used in various number theoretic computations, especially for finding primitive roots.
    pub factors: Vec<u64>,

    /// Bitmask used for fast modular reduction.
    /// Computed as 2^bit_length(modulus) - 1.
    pub mask: u64,

    /// Constants for Barrett Reduction.
    /// Used for fast modular reduction without division.
    /// b_red_constant[0] = ((2^128)/q)/(2^64)
    /// b_red_constant[1] = (2^128)/q mod 2^64
    pub b_red_constant: [u64; 2],

    /// Constant for Montgomery Reduction.
    /// Used for fast modular multiplication.
    /// m_red_constant = modulus^(-1) mod 2^64
    /// Note: This is equivalent to -modulus^(-1) mod 2^64 in 64-bit unsigned integer arithmetic.
    pub m_red_constant: u64,

    /// NTT related constants, including roots of unity and their powers.
    /// Maybe Todo: use box to store the NTTTable
    pub ntt_table: NTTTable,
}

impl SubRing {
    /// Creates a new SubRing with custom NTT transform and primitive Nth root of unity.
    ///
    /// # Arguments
    ///
    /// * `n` - Degree of the ring (must be a power of two larger than 8)
    /// * `modulus` - The modulus (should be equal to 1 modulo the root of unity)
    /// * `ntt_creator` - Function to create a custom NTT implementation
    /// * `nth_root` - The primitive Nth root of unity
    ///
    /// # Returns
    ///
    /// A Result containing either the new SubRing or an error
    pub fn new_with_custom_ntt<F>(
        n: u64,
        modulus: u64,
        ntt_creator: F,
        nth_root: u64,
    ) -> Result<Self, SubRingError>
    where
        F: FnOnce(&NTTParams) -> NTTImplementations,
    {
        // Check if N is a power of 2 and greater than the minimum
        if n < MINIMUM_RING_DEGREE_FOR_LOOP_UNROLLED_OPS || !n.is_power_of_two() {
            return Err(SubRingError::InvalidRingDegree(
                MINIMUM_RING_DEGREE_FOR_LOOP_UNROLLED_OPS,
            ));
        }

        // Check if modulus is suitable for Montgomery reduction
        // Montgomery reduction is not efficient or well-defined for certain moduli:
        // 1. When modulus is 0: Modular arithmetic is undefined.
        // 2. When modulus is a power of 2: Simpler and faster methods exist for modular arithmetic.
        //    Also, some mathematical properties required for Montgomery reduction don't hold in this case.
        // The check (modulus & (modulus - 1)) == 0 efficiently detects if modulus is a power of 2:
        // - If modulus is a power of 2, it has only one bit set.
        // - Subtracting 1 from it will set all lower bits to 1.
        // - Bitwise AND of these will be 0 only if modulus was a power of 2.
        if (modulus & (modulus - 1)) == 0 || modulus == 0 {
            return Err(SubRingError::InvalidModulusForMontgomery);
        }

        // Calculate the mask
        let mask = calculate_mask(modulus);

        // Computes Barrett reduction constants
        let b_red_constant = compute_barrett_constants(modulus);

        // Compute Montgomery reduction constant
        let m_red_constant = compute_montgomery_constant(modulus);

        let ntt = ntt_creator(&(n, modulus, nth_root, mask, b_red_constant, m_red_constant));

        Ok(Self {
            ntt,
            n: n as usize,
            modulus,
            mask,
            b_red_constant,
            m_red_constant,
            ntt_table: NTTTable::default(),
            factors: vec![],
        })
    }

    /// Evaluates p3 = p1 + p2 (mod modulus).
    #[inline(always)]
    pub fn add(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        add_vec(p1, p2, p3, self.modulus);
    }

    /// Evaluates p3 = p1 + p2.
    #[inline(always)]
    pub fn add_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        add_lazy_vec(p1, p2, p3);
    }

    /// Evaluates p3 = p1 - p2 (mod modulus).
    #[inline(always)]
    pub fn sub(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        sub_vec(p1, p2, p3, self.modulus);
    }

    /// Evaluates p3 = p1 - p2.
    #[inline(always)]
    pub fn sub_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        sub_lazy_vec(p1, p2, p3, self.modulus);
    }

    /// Evaluates p2 = -p1 (mod modulus).
    #[inline(always)]
    pub fn neg(&self, p1: &[u64], p2: &mut [u64]) {
        neg_vec(p1, p2, self.modulus);
    }

    /// Evaluates p2 = p1 (mod modulus).
    #[inline(always)]
    pub fn reduce(&self, p1: &[u64], p2: &mut [u64]) {
        reduce_vec(p1, p2, self.modulus, self.b_red_constant);
    }

    /// Evaluates p2 = p1 (mod modulus) with p2 in range [0, 2*modulus-1].
    #[inline(always)]
    pub fn reduce_lazy(&self, p1: &[u64], p2: &mut [u64]) {
        reduce_lazy_vec(p1, p2, self.modulus, self.b_red_constant);
    }

    /// Evaluates p3 = p1*p2.
    #[inline(always)]
    pub fn mul_coeffs_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_lazy_vec(p1, p2, p3);
    }

    /// Evaluates p3 = p3 + p1*p2.
    #[inline(always)]
    pub fn mul_coeffs_lazy_then_add_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_lazy_then_add_lazy_vec(p1, p2, p3);
    }

    /// Evaluates p3 = p1*p2 (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_barrett(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_barrett_vec(p1, p2, p3, self.modulus, self.b_red_constant);
    }

    /// Evaluates p3 = p1*p2 (mod modulus) with p3 in [0, 2*modulus-1].
    #[inline(always)]
    pub fn mul_coeffs_barrett_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_barrett_lazy_vec(p1, p2, p3, self.modulus, self.b_red_constant);
    }

    /// Evaluates p3 = p3 + (p1*p2) (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_barrett_then_add(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_then_add_vec(p1, p2, p3, self.modulus, self.b_red_constant);
    }

    /// Evaluates p3 = p3 + p1*p2 (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_barrett_then_add_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_barrett_then_add_lazy_vec(p1, p2, p3, self.modulus, self.b_red_constant);
    }

    /// Evaluates p3 = p1*p2 (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_montgomery(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p1*p2 (mod modulus) with p3 in range [0, 2*modulus-1].
    #[inline(always)]
    pub fn mul_coeffs_montgomery_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_lazy_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 + (p1*p2) (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_montgomery_then_add(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_then_add_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 + (p1*p2 (mod modulus)).
    #[inline(always)]
    pub fn mul_coeffs_montgomery_then_add_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_then_add_lazy_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 + p1*p2 (mod modulus) with p3 in range [0, 3modulus-2].
    #[inline(always)]
    pub fn mul_coeffs_montgomery_lazy_then_add_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_lazy_then_add_lazy_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 - p1*p2 (mod modulus).
    #[inline(always)]
    pub fn mul_coeffs_montgomery_then_sub(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_then_sub_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
    #[inline(always)]
    pub fn mul_coeffs_montgomery_then_sub_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_then_sub_lazy_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 3*modulus-2].
    #[inline(always)]
    pub fn mul_coeffs_montgomery_lazy_then_sub_lazy(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_lazy_then_sub_lazy_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
    #[inline(always)]
    pub fn mul_coeffs_montgomery_lazy_then_neg(&self, p1: &[u64], p2: &[u64], p3: &mut [u64]) {
        mul_coeffs_montgomery_lazy_then_neg_vec(p1, p2, p3, self.modulus, self.m_red_constant);
    }

    /// Evaluates p3 = (p1+p2)*scalar_mont (mod modulus).
    #[inline(always)]
    pub fn add_lazy_then_mul_scalar_montgomery(
        &self,
        p1: &[u64],
        p2: &[u64],
        scalar_mont: u64,
        p3: &mut [u64],
    ) {
        add_lazy_then_mul_scalar_montgomery_vec(
            p1,
            p2,
            scalar_mont,
            p3,
            self.modulus,
            self.m_red_constant,
        );
    }

    /// Evaluates p2 = (scalar_mont0+p1)*scalar_mont1 (mod modulus).
    #[inline(always)]
    pub fn add_scalar_lazy_then_mul_scalar_montgomery(
        &self,
        p1: &[u64],
        scalar0: u64,
        scalar_mont1: u64,
        p2: &mut [u64],
    ) {
        add_scalar_lazy_then_mul_scalar_montgomery_vec(
            p1,
            scalar0,
            scalar_mont1,
            p2,
            self.modulus,
            self.m_red_constant,
        );
    }

    /// Evaluates p2 = p1 + scalar (mod modulus).
    #[inline(always)]
    pub fn add_scalar(&self, p1: &[u64], scalar: u64, p2: &mut [u64]) {
        add_scalar_vec(p1, scalar, p2, self.modulus);
    }

    /// Evaluates p2 = p1 + scalar.
    #[inline(always)]
    pub fn add_scalar_lazy(&self, p1: &[u64], scalar: u64, p2: &mut [u64]) {
        add_scalar_lazy_vec(p1, scalar, p2);
    }

    /// Evaluates p2 = 2*modulus - p1 + scalar.
    #[inline(always)]
    pub fn add_scalar_lazy_then_neg_two_modulus_lazy(
        &self,
        p1: &[u64],
        scalar: u64,
        p2: &mut [u64],
    ) {
        add_scalar_lazy_then_neg_two_modulus_lazy_vec(p1, scalar, p2, self.modulus);
    }

    /// Evaluates p2 = p1 - scalar (mod modulus).
    #[inline(always)]
    pub fn sub_scalar(&self, p1: &[u64], scalar: u64, p2: &mut [u64]) {
        sub_scalar_vec(p1, scalar, p2, self.modulus);
    }

    /// Evaluates p2 = p1*scalar_mont (mod modulus).
    #[inline(always)]
    pub fn mul_scalar_montgomery(&self, p1: &[u64], scalar_mont: u64, p2: &mut [u64]) {
        mul_scalar_montgomery_vec(p1, scalar_mont, p2, self.modulus, self.m_red_constant);
    }

    /// Evaluates p2 = p1*scalar_mont (mod modulus) with p2 in range [0, 2*modulus-1].
    #[inline(always)]
    pub fn mul_scalar_montgomery_lazy(&self, p1: &[u64], scalar_mont: u64, p2: &mut [u64]) {
        mul_scalar_montgomery_lazy_vec(p1, scalar_mont, p2, self.modulus, self.m_red_constant);
    }

    /// Evaluates p2 = p2 + p1*scalar_mont (mod modulus).
    #[inline(always)]
    pub fn mul_scalar_montgomery_then_add(&self, p1: &[u64], scalar_mont: u64, p2: &mut [u64]) {
        mul_scalar_montgomery_then_add_vec(p1, scalar_mont, p2, self.modulus, self.m_red_constant);
    }

    /// Evaluates p2 = scalar + p1*scalar_mont (mod modulus).
    #[inline(always)]
    pub fn mul_scalar_montgomery_then_add_scalar(
        &self,
        p1: &[u64],
        scalar0: u64,
        scalar_mont1: u64,
        p2: &mut [u64],
    ) {
        mul_scalar_montgomery_then_add_scalar_vec(
            p1,
            scalar0,
            scalar_mont1,
            p2,
            self.modulus,
            self.m_red_constant,
        );
    }

    /// Evaluates p3 = (p1 + two_modulus - p2) * scalar_mont (mod modulus).
    #[inline(always)]
    pub fn sub_then_mul_scalar_montgomery_two_modulus(
        &self,
        p1: &[u64],
        p2: &[u64],
        scalar_mont: u64,
        p3: &mut [u64],
    ) {
        sub_then_mul_scalar_montgomery_two_modulus_vec(
            p1,
            p2,
            scalar_mont,
            p3,
            self.modulus,
            self.m_red_constant,
        );
    }

    /// Evaluates p2 = p1 * 2^64 (mod modulus).
    #[inline(always)]
    pub fn m_form(&self, p1: &[u64], p2: &mut [u64]) {
        m_form_vec(p1, p2, self.modulus, self.b_red_constant);
    }

    /// Evaluates p2 = p1 * 2^64 (mod modulus) with p2 in the range [0, 2*modulus-1].
    #[inline(always)]
    pub fn m_form_lazy(&self, p1: &[u64], p2: &mut [u64]) {
        m_form_lazy_vec(p1, p2, self.modulus, self.b_red_constant);
    }

    /// Evaluates p2 = p1 * (2^64)^-1 (mod modulus).    
    #[inline(always)]
    pub fn im_form(&self, p1: &[u64], p2: &mut [u64]) {
        im_form_vec(p1, p2, self.modulus, self.m_red_constant);
    }
}

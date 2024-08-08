use crate::math::{
    ntt::{
        params::{NTTParams, NTTTable},
        NTTImplementations,
    },
    ring::{
        barrett_reduction::compute_barrett_constants,
        constants::MINIMUM_RING_DEGREE_FOR_LOOP_UNROLLED_OPS,
        montgomery_reduction::compute_montgomery_constant,
    },
};

use super::errors::SubRingError;

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

    /// Creates a new SubRing with the given parameters.
    ///
    /// # Arguments
    ///
    /// * `n` - Number of coefficients (must be a power of 2)
    /// * `modulus` - The prime modulus
    ///
    /// # Returns
    ///
    /// A new SubRing instance
    pub fn new(n: usize, modulus: u64) -> Self {
        unimplemented!()
        // // Verify that n is a power of 2
        // assert!(n.is_power_of_two(), "n must be a power of 2");

        // // Compute factors of modulus - 1
        // let factors = compute_factors(modulus - 1);

        // // Compute mask
        // let mask = (1u64 << modulus.bits()) - 1;

        // // Compute Barrett reduction constants
        // let b_red_constant = compute_barrett_constants(modulus);

        // // Create NTTTable
        // let ntt_table = NTTTable::new(n, modulus);

        // // Create the NTT implementation
        // let ntt = ntt_creator();

        // SubRing {
        //     ntt,
        //     n,
        //     modulus,
        //     factors,
        //     mask,
        //     b_red_constant,
        //     m_red_constant,
        //     ntt_table,
        // }
    }

    // pub fn

    /// Performs modular addition.
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let sum = a + b;
        if sum >= self.modulus {
            sum - self.modulus
        } else {
            sum
        }
    }

    /// Performs modular subtraction.
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            self.modulus - b + a
        }
    }

    /// Performs modular multiplication using Montgomery reduction.
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.montgomery_reduce((a as u128 * b as u128) as u64)
    }

    /// Performs Montgomery reduction.
    fn montgomery_reduce(&self, t: u64) -> u64 {
        unimplemented!()

        // let m = (t as u128 * self.m_red_constant as u128) as u64;
        // let t = (t as u128 + (m as u128 * self.modulus as u128)) >> 64;
        // if t >= self.modulus {
        //     t - self.modulus
        // } else {
        //     t
        // }
    }

    // Additional methods for NTT operations, modular exponentiation, etc. would be implemented here
}

// Helper functions (implementations omitted for brevity)
fn compute_factors(n: u64) -> Vec<u64> {
    // Compute prime factors of n
    unimplemented!()
}

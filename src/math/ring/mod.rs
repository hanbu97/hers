use crate::math::ring::operations::compute_rescale_constants;

use self::{errors::RingError, subring::SubRing};
use itertools::Itertools;
use num_bigint::BigInt;

use super::ntt::NTTImplementations;

pub mod constants;
pub mod errors;
pub mod operations;
pub mod reduction;
pub mod subring;
pub mod types;

/// Ring is a structure that keeps all the variables required to operate on a polynomial represented in this ring.
#[derive(Default)]
pub struct Ring {
    /// SubRings for each level
    pub sub_rings: Vec<SubRing>,

    /// Product of the Moduli for each level
    pub modulus_at_level: Vec<BigInt>,

    /// Rescaling parameters (RNS division)
    pub rescale_constants: Vec<Vec<u64>>,

    /// Current level
    pub level: usize,
}

impl Ring {
    /// Creates a new RNS Ring with degree N and coefficient moduli Moduli with Standard NTT.
    ///
    /// # Arguments
    /// * `n` - The degree of the ring. Must be a power of two larger than 8.
    /// * `moduli` - A non-empty vector of distinct prime moduli. All moduli must also be equal to 1 modulo 2N.
    ///
    /// # Returns
    /// * `Ok(Ring)` if the parameters are valid and the Ring is successfully created.
    /// * `Err(RingError)` if the parameters are invalid.
    ///
    /// This function creates a new RNS Ring using the Standard NTT. It performs various checks on the input parameters
    /// to ensure they are valid for creating an NTT-enabling ring.
    pub fn new(n: usize, moduli: Vec<u64>) -> Result<Self, RingError> {
        unimplemented!()
    }

    pub fn new_with_custom_ntt(
        n: u64,
        moduli: Vec<u64>,
        ntt: NTTImplementations,
    ) -> Result<Self, RingError> {
        // Check if N is a power of 2
        if (n & (n - 1)) != 0 {
            return Err(RingError::InvalidRingDegree(n));
        }

        // Check if moduli is non-empty
        if moduli.is_empty() {
            return Err(RingError::EmptyModuli);
        }

        // Check if all moduli are distinct primes
        if moduli.iter().len() == moduli.iter().unique().collect::<Vec<_>>().len() {
            return Err(RingError::NonDistinctPrimeModuli);
        }

        // Compute bigQ for all levels
        let mut modulus_at_level = Vec::with_capacity(moduli.len());
        modulus_at_level.push(BigInt::from(moduli[0]));
        for &modulus in &moduli[1..] {
            let last = modulus_at_level.last().unwrap();
            modulus_at_level.push(last * BigInt::from(modulus));
        }

        // Init SubRings
        let mut sub_rings = Vec::with_capacity(moduli.len());
        for &modulus in &moduli {
            let sub_ring = SubRing::new_with_custom_ntt(n, modulus, ntt.clone())?;
            sub_rings.push(sub_ring);
        }

        // Compute rescale constants
        let rescale_constants = compute_rescale_constants(&sub_rings);

        unimplemented!()
    }
}

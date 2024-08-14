use crate::math::ring::operations::compute_rescale_constants;

use self::{errors::RingError, subring::SubRing};
use itertools::Itertools;
use num_bigint::BigInt;

use super::ntt::{
    params::{NTTParameters, NTTTable},
    NTTImplementations,
};

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
    /// * `degree` - The degree of the ring. Must be a power of two larger than 8.
    /// * `moduli` - A non-empty vector of distinct prime moduli. All moduli must also be equal to 1 modulo 2N.
    ///
    /// # Returns
    /// * `Ok(Ring)` if the parameters are valid and the Ring is successfully created.
    /// * `Err(RingError)` if the parameters are invalid.
    ///
    /// This function creates a new RNS Ring using the Standard NTT. It performs various checks on the input parameters
    /// to ensure they are valid for creating an NTT-enabling ring.
    pub fn new(degree: u64, moduli: Vec<u64>) -> anyhow::Result<Self> {
        Self::new_with_custom_ntt(
            degree,
            moduli,
            |params, table| NTTImplementations::new_standard(params, table),
            2 * degree,
        )
    }

    pub fn new_with_custom_ntt<F>(
        degree: u64,
        moduli: Vec<u64>,
        // ntt: NTTImplementations,
        ntt_creator: F,
        nth_root: u64,
    ) -> anyhow::Result<Self>
    where
        F: Fn(NTTParameters, NTTTable) -> NTTImplementations,
    {
        // Check if degree is a power of 2
        if !degree.is_power_of_two() {
            return Err(RingError::InvalidRingDegree(degree).into());
        }

        // Check if moduli is non-empty
        if moduli.is_empty() {
            return Err(RingError::EmptyModuli.into());
        }

        // Check if all moduli are distinct primes
        if moduli.iter().len() == moduli.iter().unique().collect::<Vec<_>>().len() {
            return Err(RingError::NonDistinctPrimeModuli.into());
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
            let sub_ring = SubRing::new_with_custom_ntt(degree, modulus, &ntt_creator, nth_root)?;
            sub_rings.push(sub_ring);
        }

        // Compute rescale constants
        let rescale_constants = compute_rescale_constants(&sub_rings);

        let mut ring = Ring {
            modulus_at_level,
            sub_rings,
            rescale_constants,
            level: moduli.len() - 1,
        };

        // Compute NTT constants
        ring.compute_ntt_constants(None, None)?;

        unimplemented!()
    }

    /// Computes the NTT constants for all SubRings in the Ring.
    ///
    /// This function checks that each modulus is NTT-friendly (i.e., prime and congruent to 1 mod 2N)
    /// and computes the necessary variables for the NTT.
    ///
    /// # Arguments
    ///
    /// * `primitive_roots` - Optional vector of primitive roots for each SubRing
    /// * `factors` - Optional vector of prime factors for each SubRing's modulus minus 1
    ///
    /// # Returns
    ///
    /// * `Result<(), RingError>` - Ok(()) if successful, or an error if any SubRing fails to compute its NTT constants
    pub fn compute_ntt_constants(
        &mut self,
        primitive_roots: Option<Vec<u64>>,
        factors: Option<Vec<Vec<u64>>>,
    ) -> anyhow::Result<()> {
        for (i, sub_ring) in self.sub_rings.iter_mut().enumerate() {
            if let (Some(roots), Some(facts)) = (&primitive_roots, &factors) {
                sub_ring.ntt_table.primitive_root = roots[i];
                sub_ring.factors = facts[i].clone();
            }

            sub_ring.compute_ntt_constants()?;
        }

        Ok(())
    }
}

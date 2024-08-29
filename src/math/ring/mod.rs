use crate::math::ring::operations::compute_rescale_constants;

use self::{errors::RingError, polynomial::Poly, subring::SubRing};
use itertools::Itertools;
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

use super::ntt::{
    params::{NTTParameters, NTTTable},
    NTTImplementations,
};

pub mod constants;
pub mod errors;
pub mod operations;
pub mod polynomial;
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
            |params, _table| NTTImplementations::new_standard(params),
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
        if moduli.iter().unique().count() < moduli.len() {
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

        Ok(ring)
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

    /// NewPoly creates a new polynomial with all coefficients set to 0.
    pub fn new_poly(&self) -> Poly {
        Poly::new(self.degree() as usize, self.level)
    }

    // Operations
    /// Add evaluates p3 = p1 + p2 coefficient-wise in the ring.
    pub fn add(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].add(&p1.coeffs[i], &p2.coeffs[i], &mut p3.coeffs[i]);
        }
    }

    /// AddLazy evaluates p3 = p1 + p2 coefficient-wise in the ring, with p3 in [0, 2*modulus-1].
    pub fn add_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].add_lazy(&p1.coeffs[i], &p2.coeffs[i], &mut p3.coeffs[i]);
        }
    }

    /// Sub evaluates p3 = p1 - p2 coefficient-wise in the ring.
    pub fn sub(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].sub(&p1.coeffs[i], &p2.coeffs[i], &mut p3.coeffs[i]);
        }
    }

    /// Neg evaluates p2 = -p1 coefficient-wise in the ring.
    pub fn neg(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].neg(&p1.coeffs[i], &mut p2.coeffs[i]);
        }
    }

    /// Reduce evaluates p2 = p1 coefficient-wise mod modulus in the ring.
    pub fn reduce(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].reduce(&p1.coeffs[i], &mut p2.coeffs[i]);
        }
    }

    /// ReduceLazy evaluates p2 = p1 coefficient-wise mod modulus in the ring, with p2 in [0, 2*modulus-1].
    pub fn reduce_lazy(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].reduce_lazy(&p1.coeffs[i], &mut p2.coeffs[i]);
        }
    }

    /// SubLazy evaluates p3 = p1 - p2 coefficient-wise in the ring, with p3 in [0, 2*modulus-1].
    pub fn sub_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].sub_lazy(&p1.coeffs[i], &p2.coeffs[i], &mut p3.coeffs[i]);
        }
    }

    /// MulCoeffsBarrett evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction.
    pub fn mul_coeffs_barrett(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_barrett(&p1.coeffs[i], &p2.coeffs[i], &mut p3.coeffs[i]);
        }
    }

    /// MulCoeffsBarrettLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_barrett_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_barrett_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsBarrettThenAdd evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Barrett reduction.
    pub fn mul_coeffs_barrett_then_add(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_barrett_then_add(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsBarrettThenAddLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_barrett_then_add_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_barrett_then_add_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomery evaluates p3 = p1 * p2 coefficient-wise in the ring, with Montgomery reduction.
    pub fn mul_coeffs_montgomery(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryLazyThenNeg evaluates p3 = -p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_lazy_then_neg(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_lazy_then_neg(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryThenAdd evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_then_add(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_then_add(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryThenAddLazy evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_then_add_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_then_add_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryThenSub evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction.
    pub fn mul_coeffs_montgomery_then_sub(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_then_sub(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryThenSubLazy evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_then_sub_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_then_sub_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    /// MulCoeffsMontgomeryThenSubLazy evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
    pub fn mul_coeffs_montgomery_lazy_then_sub_lazy(&self, p1: &Poly, p2: &Poly, p3: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].mul_coeffs_montgomery_lazy_then_sub_lazy(
                &p1.coeffs[i],
                &p2.coeffs[i],
                &mut p3.coeffs[i],
            );
        }
    }

    // AddScalar evaluates p2 = p1 + scalar coefficient-wise in the ring.
    pub fn add_scalar(&self, p1: &Poly, scalar: u64, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].add_scalar(&p1.coeffs[i], scalar, &mut p2.coeffs[i]);
        }
    }

    // AddScalarBigint evaluates p2 = p1 + scalar coefficient-wise in the ring.
    pub fn add_scalar_bigint(&self, p1: &Poly, scalar: &BigInt, p2: &mut Poly) {
        for i in 0..=self.level {
            let tmp = scalar % BigInt::from(self.sub_rings[i].modulus);
            self.sub_rings[i].add_scalar(&p1.coeffs[i], tmp.to_u64().unwrap(), &mut p2.coeffs[i]);
        }
    }

    /// Returns the list of primes in the modulus chain.
    pub fn moduli_chain(&self) -> Vec<u64> {
        self.sub_rings.iter().map(|sr| sr.modulus).collect()
    }

    /// ModuliChainLength returns the number of primes in the RNS basis of the ring.
    pub fn moduli_chain_length(&self) -> usize {
        self.sub_rings.len()
    }

    /// Returns the ring degree.
    pub fn degree(&self) -> u64 {
        self.sub_rings[0].degree
    }

    /// Evaluates p2 = NTT(p1)
    pub fn ntt(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].ntt(
                &p1.coeffs[i],
                &mut p2.coeffs[i],
                &self.sub_rings[i].ntt_table,
            );
        }
    }

    /// Evaluates p2 = p1 * (2^64)^-1 coefficient-wise in the ring.
    pub fn m_form(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].m_form(&p1.coeffs[i], &mut p2.coeffs[i]);
        }
    }
}

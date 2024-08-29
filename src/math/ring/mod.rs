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

    /// INTT evaluates p2 = INTT(p1).
    pub fn intt(&self, p1: &Poly, p2: &mut Poly) {
        for i in 0..=self.level {
            self.sub_rings[i].intt(
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

    // Shift evaluates p2 = p2<<<k coefficient-wise in the ring.
    pub fn shift(&self, p1: &Poly, k: usize, p2: &mut Poly) {
        for (src, dst) in p1.coeffs.iter().zip(p2.coeffs.iter_mut()) {
            dst.copy_from_slice(src);
            dst.rotate_left(k);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_new_ring() -> Result<()> {
        // Test with invalid degree (not power of 2)
        assert!(Ring::new(3, vec![97]).is_err());

        // Test with empty moduli
        assert!(Ring::new(16, vec![]).is_err());

        // Test with non-distinct moduli
        assert!(Ring::new(16, vec![97, 97]).is_err());

        // Test with non-prime modulus
        let ring = Ring::new(16, vec![4]);
        assert!(ring.is_err());

        // Test with non NTT-friendly modulus
        let ring = Ring::new(16, vec![7]);
        assert!(ring.is_err());

        // Test with one NTT-friendly and one non NTT-friendly modulus
        let ring = Ring::new(16, vec![97, 7]);
        assert!(ring.is_err());

        // Test with valid parameters
        let ring = Ring::new(16, vec![97]);
        assert!(ring.is_ok());

        Ok(())
    }

    #[test]
    fn test_shift() {
        let r = Ring::new(16, vec![97]).unwrap();
        let mut p1 = r.new_poly();
        let mut p2 = r.new_poly();

        // Initialize p1 with values 0 to 15
        for i in 0..16 {
            p1.coeffs[0][i] = i as u64;
        }

        r.shift(&p1, 3, &mut p2);

        assert_eq!(
            p2.coeffs[0],
            vec![3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2]
        );
    }

    struct Parameters {
        pub log_n: u64,
        pub qi: Vec<u64>,
        pub pi: Vec<u64>,
    }

    const T: u64 = 0x3ee0001;
    const DEFAULT_SIGMA: f64 = 3.2;
    const DEFAULT_BOUND: f64 = 6.0 * DEFAULT_SIGMA;

    // struct TestParams {
    //     ring_q: Ring,
    //     ring_p: Ring,
    //     prng: PRNG,
    //     uniform_sampler_q: UniformSampler,
    //     uniform_sampler_p: UniformSampler,
    // }

    const QI60: [u64; 32] = [
        0x1fffffffffe00001,
        0x1fffffffffc80001,
        0x1fffffffffb40001,
        0x1fffffffff500001,
        0x1fffffffff380001,
        0x1fffffffff000001,
        0x1ffffffffef00001,
        0x1ffffffffee80001,
        0x1ffffffffeb40001,
        0x1ffffffffe780001,
        0x1ffffffffe600001,
        0x1ffffffffe4c0001,
        0x1ffffffffdf40001,
        0x1ffffffffdac0001,
        0x1ffffffffda40001,
        0x1ffffffffc680001,
        0x1ffffffffc000001,
        0x1ffffffffb880001,
        0x1ffffffffb7c0001,
        0x1ffffffffb300001,
        0x1ffffffffb1c0001,
        0x1ffffffffadc0001,
        0x1ffffffffa400001,
        0x1ffffffffa140001,
        0x1ffffffff9d80001,
        0x1ffffffff9140001,
        0x1ffffffff8ac0001,
        0x1ffffffff8a80001,
        0x1ffffffff81c0001,
        0x1ffffffff7800001,
        0x1ffffffff7680001,
        0x1ffffffff7080001,
    ];

    const PI60: [u64; 32] = [
        0x1ffffffff6c80001,
        0x1ffffffff6140001,
        0x1ffffffff5f40001,
        0x1ffffffff5700001,
        0x1ffffffff4bc0001,
        0x1ffffffff4380001,
        0x1ffffffff3240001,
        0x1ffffffff2dc0001,
        0x1ffffffff1a40001,
        0x1ffffffff11c0001,
        0x1ffffffff0fc0001,
        0x1ffffffff0d80001,
        0x1ffffffff0c80001,
        0x1ffffffff08c0001,
        0x1fffffffefd00001,
        0x1fffffffef9c0001,
        0x1fffffffef600001,
        0x1fffffffeef40001,
        0x1fffffffeed40001,
        0x1fffffffeed00001,
        0x1fffffffeebc0001,
        0x1fffffffed540001,
        0x1fffffffed440001,
        0x1fffffffed2c0001,
        0x1fffffffed200001,
        0x1fffffffec940001,
        0x1fffffffec6c0001,
        0x1fffffffebe80001,
        0x1fffffffebac0001,
        0x1fffffffeba40001,
        0x1fffffffeb4c0001,
        0x1fffffffeb280001,
    ];

    lazy_static::lazy_static! {
        static ref TEST_PARAMETERS: Vec<Parameters> = vec![Parameters {
            log_n: 10,
            qi: QI60[QI60.len() - 14..].to_vec(),
            pi: PI60[PI60.len() - 14..].to_vec(),
        },];
    }

    fn test_string(opname: &str, ring_q: &Ring) -> String {
        format!(
            "{}/N={}/limbs={}",
            opname,
            ring_q.degree(),
            ring_q.moduli_chain_length()
        )
    }
}

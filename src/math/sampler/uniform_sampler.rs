use rand_core::RngCore;

use crate::math::ring::{polynomial::Poly, Ring};

pub struct UniformSampler<R: RngCore + Clone> {
    base_ring: Ring,
    prng: R,
}

impl<R: RngCore + Clone> UniformSampler<R> {
    pub fn new(prng: R, base_ring: Ring) -> Self {
        UniformSampler { base_ring, prng }
    }

    pub fn at_level(&self, level: usize) -> Self {
        UniformSampler {
            base_ring: self.base_ring.at_level(level),
            prng: self.prng.clone(),
        }
    }

    pub fn read(&mut self, pol: &mut Poly) {
        self.read_inner(pol, |_, b, _| b);
    }

    pub fn read_and_add(&mut self, pol: &mut Poly) {
        self.read_inner(pol, |a, b, c| (a + b) % c);
    }

    fn read_inner<F>(&mut self, pol: &mut Poly, f: F)
    where
        F: Fn(u64, u64, u64) -> u64,
    {
        let level = self.base_ring.moduli_chain().len() - 1;
        let n = self.base_ring.degree();

        for j in 0..=level {
            let qi = self.base_ring.moduli_chain()[j];
            let mask = (1u64 << (qi.ilog2() + 1)) - 1;

            for i in 0..n {
                let mut random_uint = self.prng.next_u64() & mask;
                while random_uint >= qi {
                    random_uint = self.prng.next_u64() & mask;
                }

                pol.coeffs[j][i as usize] = f(pol.coeffs[j][i as usize], random_uint, qi);
            }
        }
    }

    pub fn read_new(&mut self) -> Poly {
        let mut pol = self.base_ring.new_poly();
        self.read(&mut pol);
        pol
    }
}

pub fn rand_uniform<R: RngCore>(prng: &mut R, v: u64, mask: u64) -> u64 {
    loop {
        let random_int = rand_int64(prng, mask);
        if random_int < v {
            return random_int;
        }
    }
}

fn rand_int64<R: RngCore>(prng: &mut R, mask: u64) -> u64 {
    let mut random_bytes = [0u8; 8];
    prng.fill_bytes(&mut random_bytes);
    u64::from_be_bytes(random_bytes) & mask
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::rand::fixed_rand::FixedRandom;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_uniform_sampler() {
        let ring = Ring::new(16, vec![97]).unwrap();
        let prng = ChaCha20Rng::seed_from_u64(12345);
        let mut sampler = UniformSampler::new(prng, ring);

        let pol = sampler.read_new();

        // Basic sanity checks
        assert_eq!(pol.coeffs.len(), 1);
        assert_eq!(pol.coeffs[0].len(), 16);

        // Check if coefficients are within bounds
        for coeff in &pol.coeffs[0] {
            assert!(*coeff < 97);
        }
    }

    #[test]
    fn test_uniform_sampler_for_go_comparison() {
        let degree = 16;
        let modulus = 97u64;

        let ring = Ring::new(degree, vec![modulus]).unwrap();

        let random_numbers: Vec<u64> = vec![
            12300691266768169957,
            11820962631477129540,
            17883164048736053457,
            3763972486863574293,
            1820076341489734165,
            3641032664938027597,
            13856849471837345700,
            11093809191189476932,
            3825727900517507516,
            1833057533073052521,
            1862332976669812221,
            839969800260750646,
            5102443919682075139,
            3469586806049992227,
            1289589544690565689,
            4686543058686886181,
            4863912002469269231,
            5872498502673663872,
            10304755672215287095,
            1191274978463715656,
        ];
        let fixed_rng = FixedRandom::new(random_numbers.clone());

        let mut sampler = UniformSampler::new(fixed_rng, ring.clone());

        let mut pol = ring.new_poly();

        sampler.read(&mut pol);

        let coeffs = &pol.coeffs[0];

        let expected_coeffs = vec![68, 81, 21, 21, 77, 36, 68, 60, 54, 3, 35, 57, 37, 0, 55, 72];
        assert_eq!(
            coeffs, &expected_coeffs,
            "Sampled coefficients do not match Go implementation"
        );

        for (i, &coeff) in coeffs.iter().enumerate() {
            assert!(
                coeff < modulus,
                "Coefficient at index {} is out of range: {}",
                i,
                coeff
            );
        }
    }
}

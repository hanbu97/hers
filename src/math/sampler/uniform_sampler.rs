use std::cell::Cell;

use crate::math::ring::{polynomial::Poly, reduction::conditional::c_red, Ring};

use super::{traits::Sampler, KeyedPRNG};

/// Uniform sampler for polynomials.
pub struct UniformSampler {
    base_ring: Ring,
    prng: KeyedPRNG,
    random_buffer: Vec<u8>,
    ptr: Cell<usize>,
}

impl UniformSampler {
    /// Creates a new UniformSampler.
    ///
    /// # Arguments
    ///
    /// * `prng` - A keyed pseudo-random number generator.
    /// * `base_ring` - The ring in which polynomials will be sampled.
    pub fn new(prng: KeyedPRNG, base_ring: Ring) -> Self {
        let n = base_ring.degree() as usize;
        UniformSampler {
            base_ring,
            prng,
            random_buffer: vec![0; n.max(1024)],
            ptr: Cell::new(0),
        }
    }

    /// Internal method for reading or adding sampled values.
    ///
    /// # Arguments
    ///
    /// * `pol` - The polynomial to be filled or added to.
    /// * `f` - A function that determines how to combine existing and new values.
    fn read_internal<F>(&mut self, pol: &mut Poly, f: F)
    where
        F: Fn(u64, u64, u64) -> u64,
    {
        let level = self.base_ring.level;
        let n = self.base_ring.degree() as usize;

        let mut ptr = self.ptr.get();

        for j in 0..=level {
            let qi = self.base_ring.sub_rings[j].modulus;
            let mask = self.base_ring.sub_rings[j].mask;

            let coeffs = &mut pol.coeffs[j];

            for i in 0..n {
                let random_uint = loop {
                    if ptr == 0 || ptr == n {
                        self.prng.read(&mut self.random_buffer).unwrap();
                        ptr = 0;
                    }

                    let random_uint =
                        u64::from_be_bytes(self.random_buffer[ptr..ptr + 8].try_into().unwrap())
                            & mask;
                    ptr += 8;

                    if random_uint < qi {
                        break random_uint;
                    }
                };

                coeffs[i] = f(coeffs[i], random_uint, qi);
            }
        }

        self.ptr.set(ptr);
    }

    /// Creates a new UniformSampler with a different PRNG.
    pub fn with_prng(&self, prng: KeyedPRNG) -> Self {
        UniformSampler {
            base_ring: self.base_ring.clone(),
            prng,
            random_buffer: self.random_buffer.clone(),
            ptr: Cell::new(0),
        }
    }
}

impl Sampler for UniformSampler {
    fn read(&mut self, pol: &mut Poly) {
        self.read_internal(pol, |_, b, _| b);
    }

    fn read_new(&mut self) -> Poly {
        let mut pol = self.base_ring.new_poly();
        self.read(&mut pol);
        pol
    }

    fn read_and_add(&mut self, pol: &mut Poly) {
        self.read_internal(pol, |a, b, c| c_red(a.wrapping_add(b), c));
    }

    fn at_level(&self, level: usize) -> Self {
        UniformSampler {
            base_ring: self.base_ring.at_level(level),
            prng: KeyedPRNG::new(&self.prng.key()),
            random_buffer: self.random_buffer.clone(),
            ptr: Cell::new(0),
        }
    }
}

/// Samples a uniform random integer in the range [0, v-1].
pub fn rand_uniform(prng: &mut KeyedPRNG, v: u64, mask: u64) -> u64 {
    loop {
        let random_int = rand_int64(prng, mask);
        if random_int < v {
            return random_int;
        }
    }
}

/// Samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 32].
pub fn rand_int32(prng: &mut KeyedPRNG, mask: u64) -> u64 {
    let mut random_bytes = [0u8; 4];
    prng.read(&mut random_bytes)
        .expect("Failed to read from PRNG");
    let random_uint32 = u32::from_be_bytes(random_bytes) as u64;
    mask & random_uint32
}

/// Samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
pub fn rand_int64(prng: &mut KeyedPRNG, mask: u64) -> u64 {
    let mut random_bytes = [0u8; 8];
    prng.read(&mut random_bytes)
        .expect("Failed to read from PRNG");
    let random_uint64 = u64::from_be_bytes(random_bytes);
    mask & random_uint64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rand_uniform() {
        let mut prng = KeyedPRNG::new(&[0u8; 32]);
        let v = 100;
        let mask = 0x7f; // 2^7 - 1

        for _ in 0..1000 {
            let random_int = rand_uniform(&mut prng, v, mask);
            assert!(random_int < v);
        }
    }

    #[test]
    fn test_rand_int32() {
        let mut prng = KeyedPRNG::new(&[0u8; 32]);
        let mask = 0xffffffff; // 2^32 - 1

        for _ in 0..1000 {
            let random_int = rand_int32(&mut prng, mask);
            assert!(random_int <= mask);
        }
    }

    #[test]
    fn test_rand_int64() {
        let mut prng = KeyedPRNG::new(&[0u8; 32]);
        let mask = 0xffffffffffffffff; // 2^64 - 1

        for _ in 0..1000 {
            let random_int = rand_int64(&mut prng, mask);
            assert!(random_int <= mask);
        }
    }
}

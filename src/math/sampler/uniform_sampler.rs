// use super::{
//     traits::{Ring, Sampler},
//     KeyedPRNG,
// };
// use std::cell::Cell;

// pub struct UniformSampler<R: Ring> {
//     base_ring: R,
//     prng: KeyedPRNG,
//     random_buffer: Vec<u8>,
//     ptr: Cell<usize>,
// }

// impl<R: Ring> UniformSampler<R> {
//     pub fn new(prng: KeyedPRNG, base_ring: R) -> Self {
//         let n = base_ring.n();
//         UniformSampler {
//             base_ring,
//             prng,
//             random_buffer: vec![0; n.max(1024)],
//             ptr: Cell::new(0),
//         }
//     }

//     fn read_internal<F>(&mut self, pol: &mut R::Poly, f: F)
//     where
//         F: Fn(u64, u64, u64) -> u64,
//     {
//         let level = self.base_ring.level();
//         let n = self.base_ring.n();

//         let mut ptr = self.ptr.get();
//         if ptr == 0 || ptr == n {
//             self.prng.read(&mut self.random_buffer).unwrap();
//             ptr = 0;
//         }

//         let buffer = &self.random_buffer;

//         for j in 0..=level {
//             let qi = self.base_ring.modulus(j);
//             let mask = self.base_ring.mask(j);

//             let coeffs = pol.coeffs_mut(j);

//             for i in 0..n {
//                 let mut random_uint = loop {
//                     if ptr == n {
//                         self.prng.read(&mut self.random_buffer).unwrap();
//                         ptr = 0;
//                     }

//                     let random_uint =
//                         u64::from_be_bytes(buffer[ptr..ptr + 8].try_into().unwrap()) & mask;
//                     ptr += 8;

//                     if random_uint < qi {
//                         break random_uint;
//                     }
//                 };

//                 coeffs[i] = f(coeffs[i], random_uint, qi);
//             }
//         }

//         self.ptr.set(ptr);
//     }
// }

// impl<R: Ring> Sampler for UniformSampler<R> {
//     type Ring = R;

//     fn read(&mut self, pol: &mut R::Poly) {
//         self.read_internal(pol, |_, b, _| b);
//     }

//     fn read_new(&mut self) -> R::Poly {
//         let mut pol = self.base_ring.new_poly();
//         self.read(&mut pol);
//         pol
//     }

//     fn read_and_add(&mut self, pol: &mut R::Poly) {
//         self.read_internal(pol, |a, b, c| cred(a.wrapping_add(b), c));
//     }

//     fn at_level(&self, level: usize) -> Self {
//         UniformSampler {
//             base_ring: self.base_ring.at_level(level),
//             prng: KeyedPRNG::new(&self.prng.key()),
//             random_buffer: self.random_buffer.clone(),
//             ptr: Cell::new(0),
//         }
//     }
// }

// // Helper function for Chinese Remainder reduction
// fn cred(a: u64, m: u64) -> u64 {
//     if a >= m {
//         a % m
//     } else {
//         a
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::traits::Ring;

//     struct TestRing {
//         n: usize,
//         moduli: Vec<u64>,
//         level: usize,
//     }

//     impl Ring for TestRing {
//         type Poly = Vec<Vec<u64>>;

//         fn new(n: usize, moduli: Vec<u64>) -> Self {
//             TestRing {
//                 n,
//                 moduli: moduli.clone(),
//                 level: moduli.len() - 1,
//             }
//         }

//         fn n(&self) -> usize {
//             self.n
//         }

//         fn modulus(&self, i: usize) -> u64 {
//             self.moduli[i]
//         }

//         fn level(&self) -> usize {
//             self.level
//         }

//         fn mask(&self, i: usize) -> u64 {
//             (1 << 64) - 1
//         }

//         fn new_poly(&self) -> Self::Poly {
//             vec![vec![0; self.n]; self.level + 1]
//         }

//         fn at_level(&self, level: usize) -> Self {
//             TestRing {
//                 n: self.n,
//                 moduli: self.moduli[..=level].to_vec(),
//                 level,
//             }
//         }
//     }

//     #[test]
//     fn test_uniform_sampler() {
//         let ring = TestRing::new(1024, vec![0xffffff01, 0xffffc001]);
//         let prng = KeyedPRNG::new_random();
//         let mut sampler = UniformSampler::new(prng, ring);

//         let mut pol = sampler.read_new();

//         // Check that all coefficients are within the correct range
//         for level in 0..=sampler.base_ring.level() {
//             let modulus = sampler.base_ring.modulus(level);
//             for coeff in &pol[level] {
//                 assert!(*coeff < modulus);
//             }
//         }

//         // Test read_and_add
//         let original = pol.clone();
//         sampler.read_and_add(&mut pol);

//         // Check that coefficients have changed and are still within range
//         for level in 0..=sampler.base_ring.level() {
//             let modulus = sampler.base_ring.modulus(level);
//             for (orig, new) in original[level].iter().zip(pol[level].iter()) {
//                 assert_ne!(orig, new);
//                 assert!(*new < modulus);
//             }
//         }
//     }
// }

use blake3::Hasher;
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

pub mod traits;

pub mod gaussian_sampler;
pub mod sampler;
pub mod ternary_sampler;
pub mod uniform_sampler;

use std::io;

pub trait PRNG: io::Read {}

impl<T: io::Read> PRNG for T {}

/// KeyedPRNG is a structure storing the parameters used to securely and deterministically generate shared
/// sequences of random bytes among different parties using the hash function blake3. Backward sequence
/// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
/// security (given the digest i, compute the digest i+1) is only ensured if the KeyedPRNG is keyed.
pub struct KeyedPRNG {
    key: [u8; 32],
    hasher: Hasher,
}

impl KeyedPRNG {
    /// NewKeyedPRNG creates a new instance of KeyedPRNG.
    /// Accepts a key, which is padded or truncated to 32 bytes.
    /// WARNING: A PRNG INITIALISED WITH an empty key IS INSECURE!
    pub fn new(key: &[u8]) -> Self {
        let mut padded_key = [0u8; 32];
        let len = std::cmp::min(key.len(), 32);
        padded_key[..len].copy_from_slice(&key[..len]);

        let hasher = Hasher::new_keyed(&padded_key);
        KeyedPRNG {
            key: padded_key,
            hasher,
        }
    }

    /// NewPRNG creates KeyedPRNG keyed from a cryptographically secure random number generator
    /// for instances where no key should be provided by the user
    pub fn new_random() -> Self {
        let mut rng = ChaCha20Rng::from_entropy();
        let mut key = [0u8; 32];
        rng.fill_bytes(&mut key);
        Self::new(&key)
    }

    /// Reset resets the PRNG to its initial state.
    pub fn reset(&mut self) {
        self.hasher = Hasher::new_keyed(&self.key);
    }

    /// Key returns a copy of the key used to seed the PRNG.
    /// This value can be used with `new` to instantiate
    /// a new PRNG that will produce the same stream of bytes.
    pub fn key(&self) -> [u8; 32] {
        self.key
    }

    /// Read reads bytes from the KeyedPRNG on sum.
    pub fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut xof = self.hasher.finalize_xof();
        xof.fill(buf);
        Ok(buf.len())
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_prng() {
        let key = [
            0x49, 0x0a, 0x42, 0x3d, 0x97, 0x9d, 0xc1, 0x07, 0xa1, 0xd7, 0xe9, 0x7b, 0x3b, 0xce,
            0xa1, 0xdb, 0x42, 0xf3, 0xa6, 0xd5, 0x75, 0xd2, 0x0c, 0x92, 0xb7, 0x35, 0xce, 0x0c,
            0xee, 0x09, 0x7c, 0x98,
        ];

        let mut prng_a = KeyedPRNG::new(&key);
        let mut prng_b = KeyedPRNG::new(&key);

        let mut sum0 = vec![0u8; 512];
        let mut sum1 = vec![0u8; 512];

        prng_b.read(&mut sum1).unwrap();

        for _ in 0..128 {
            prng_b.read(&mut sum1).unwrap();
        }

        prng_b.reset();

        prng_a.read(&mut sum0).unwrap();
        prng_b.read(&mut sum1).unwrap();

        assert_eq!(sum0, sum1);
    }

    #[test]
    fn test_prng_with_short_key() {
        let short_key = [0x49, 0x0a, 0x42, 0x3d];
        let mut prng = KeyedPRNG::new(&short_key);
        let mut output = vec![0u8; 32];
        prng.read(&mut output).unwrap();
        assert!(!output.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_prng_random() {
        let mut prng = KeyedPRNG::new_random();
        let mut output1 = vec![0u8; 32];
        let mut output2 = vec![0u8; 32];
        prng.read(&mut output1).unwrap();
        prng.reset();
        prng.read(&mut output2).unwrap();
        assert_eq!(output1, output2);
    }
}

use blake3::Hasher;
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

pub struct KeyedPRNG {
    hasher: Hasher,
}

impl KeyedPRNG {
    pub fn new(key: &[u8]) -> Self {
        let mut padded_key = [0u8; 32];
        let len = std::cmp::min(key.len(), 32);
        padded_key[..len].copy_from_slice(&key[..len]);

        let hasher = Hasher::new_keyed(&padded_key);
        KeyedPRNG { hasher }
    }

    pub fn new_random() -> Self {
        let mut rng = ChaCha20Rng::from_entropy();
        let mut key = [0u8; 32];
        rng.fill_bytes(&mut key);
        Self::new(&key)
    }

    pub fn read(&mut self, buf: &mut [u8]) {
        let mut xof = self.hasher.finalize_xof();
        xof.fill(buf);
    }

    pub fn reset(&mut self) {
        self.hasher.reset();
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

        prng_b.read(&mut sum1);

        for _ in 0..128 {
            prng_b.read(&mut sum1);
        }

        prng_b.reset();

        prng_a.read(&mut sum0);
        prng_b.read(&mut sum1);

        assert_eq!(sum0, sum1);
    }

    #[test]
    fn test_prng_with_short_key() {
        let short_key = [0x49, 0x0a, 0x42, 0x3d];
        let mut prng = KeyedPRNG::new(&short_key);
        let mut output = vec![0u8; 32];
        prng.read(&mut output);
        assert!(!output.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_prng_random() {
        let mut prng = KeyedPRNG::new_random();
        let mut output1 = vec![0u8; 32];
        let mut output2 = vec![0u8; 32];
        prng.read(&mut output1);
        prng.reset();
        prng.read(&mut output2);
        assert_eq!(output1, output2);
    }
}

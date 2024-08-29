// use rug::rand::RandState;
// use rug::Assign;
// use rug::Integer;

// pub fn rand_u64(rng: &mut RandState) -> u64 {
//     rng.bits(64) as u64
// }

// pub fn rand_f64(rng: &mut RandState, min: f64, max: f64) -> f64 {
//     min + (max - min) * (rng.bits(53) as f64 / (1u64 << 53) as f64)
// }

// pub fn rand_complex(rng: &mut RandState, min: f64, max: f64) -> (f64, f64) {
//     (rand_f64(rng, min, max), rand_f64(rng, min, max))
// }

// pub fn rand_int(rng: &mut RandState, max: &Integer) -> Integer {
//     let mut result = Integer::new();
//     result.assign(max);
//     result.random_below(rng)
// }

// pub struct KeyedPRNG<'a> {
//     rand_state: RandState<'a>,
// }

// impl<'a> KeyedPRNG<'a> {
//     // pub fn new(key: Option<&Integer>) -> Self {
//     //     let mut rand_state = RandState::new();
//     //     if let Some(k) = key {
//     //         rand_state.seed(k);
//     //     }
//     //     KeyedPRNG { rand_state }
//     // }

//     pub fn new(key: &[u8]) -> Self {
//         let mut rand_state = RandState::new();
//         let seed = Integer::from_digits(key, rug::integer::Order::Lsf);
//         rand_state.seed(&seed);
//         KeyedPRNG { rand_state }
//     }

//     pub fn read(&mut self, buf: &mut [u8]) {
//         for chunk in buf.chunks_mut(4) {
//             let random_bits = self.rand_state.bits(chunk.len() as u32 * 8);
//             for (i, byte) in chunk.iter_mut().enumerate() {
//                 *byte = (random_bits >> (i * 8)) as u8;
//             }
//         }
//     }

//     pub fn reset(&mut self) {
//         self.rand_state = RandState::new();
//     }
// }

// impl std::io::Read for KeyedPRNG<'_> {
//     fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
//         self.read(buf);
//         Ok(buf.len())
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_prng() {
//         let key: [u8; 32] = [
//             0x49, 0x0a, 0x42, 0x3d, 0x97, 0x9d, 0xc1, 0x07, 0xa1, 0xd7, 0xe9, 0x7b, 0x3b, 0xce,
//             0xa1, 0xdb, 0x42, 0xf3, 0xa6, 0xd5, 0x75, 0xd2, 0x0c, 0x92, 0xb7, 0x35, 0xce, 0x0c,
//             0xee, 0x09, 0x7c, 0x98,
//         ];

//         let mut ha = KeyedPRNG::new(&key);
//         let mut hb = KeyedPRNG::new(&key);

//         let mut sum0 = vec![0u8; 512];
//         let mut sum1 = vec![0u8; 512];

//         for _ in 0..128 {
//             hb.read(&mut sum1);
//         }

//         hb.reset();

//         ha.read(&mut sum0);
//         hb.read(&mut sum1);

//         assert_eq!(sum0, sum1);
//     }
// }

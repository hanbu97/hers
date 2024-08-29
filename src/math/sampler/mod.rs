#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prng() {
        let key: [u8; 32] = [
            0x49, 0x0a, 0x42, 0x3d, 0x97, 0x9d, 0xc1, 0x07, 0xa1, 0xd7, 0xe9, 0x7b, 0x3b, 0xce,
            0xa1, 0xdb, 0x42, 0xf3, 0xa6, 0xd5, 0x75, 0xd2, 0x0c, 0x92, 0xb7, 0x35, 0xce, 0x0c,
            0xee, 0x09, 0x7c, 0x98,
        ];
    }
}

use blake2::{
    digest::{Update, VariableOutput},
    Blake2bVarCore,
};

pub struct KeyedPRNG {
    xof: Blake2bVarCore,
}

// impl KeyedPRNG {
//     pub fn new(key: &[u8]) -> anyhow::Result<Self> {
//         let xof = Blake2bVarCore::new_with_params(&[], &[], key.len(), 64);
//         let mut prng = KeyedPRNG { xof };
//         prng.xof.update(key);
//         Ok(prng)
//     }

//     pub fn read(&mut self, sum: &mut [u8]) -> Result<(), blake2::Error> {
//         self.xof.finalize_xof_into(sum);
//         Ok(())
//     }

//     pub fn reset(&mut self) {
//         // Reset the internal state
//         self.xof = Blake2bVarCore::new_with_params(&[], &[], 0, 64).unwrap();
//     }
// }

// Example usage
// fn main() -> Result<(), blake2::Error> {
//     let key = [0x49, 0x0a, 0x42, 0x3d, 0x97, 0x9d, 0xc1, 0x07];
//     let mut prng = KeyedPRNG::new(&key)?;

//     let mut sum = [0u8; 512];
//     prng.read(&mut sum)?;
//     println!("Random bytes: {:?}", &sum[..10]); // Print first 10 bytes

//     prng.reset();

//     let mut sum2 = [0u8; 512];
//     prng.read(&mut sum2)?;
//     println!("Random bytes after reset: {:?}", &sum2[..10]);

//     assert_eq!(sum, sum2, "Outputs should be the same after reset");

//     Ok(())
// }

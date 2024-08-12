use super::*;

pub trait BigUintExt {
    fn get_limb(&self, i: u64) -> u64;
}

impl BigUintExt for BigUint {
    fn get_limb(&self, i: u64) -> u64 {
        let digits = self.to_u64_digits();
        if i >= digits.len() as u64 {
            0
        } else {
            digits[i as usize]
        }
    }
}

impl BigUintExt for BigInt {
    fn get_limb(&self, i: u64) -> u64 {
        // We only care about the magnitude for limbs, not the sign
        let digits = self.magnitude().to_u64_digits();
        if i >= digits.len() as u64 {
            0
        } else {
            digits[i as usize]
        }
    }
}

// pub trait BigUintExt {
//     fn get_limb(&self, i: u64) -> u32;
// }

// impl BigUintExt for BigUint {
//     fn get_limb(&self, i: u64) -> u32 {
//         let digits = self.to_u32_digits();
//         if i >= digits.len() as u64 {
//             0
//         } else {
//             digits[i as usize]
//         }
//     }
// }

// impl BigUintExt for BigInt {
//     fn get_limb(&self, i: u64) -> u32 {
//         // We only care about the magnitude for limbs, not the sign
//         let digits = self.magnitude().to_u32_digits();
//         if i >= digits.len() as u64 {
//             0
//         } else {
//             digits[i as usize]
//         }
//     }
// }

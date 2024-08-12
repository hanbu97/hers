use self::constants::BITS;
use self::ext::BigUintExt;

use super::*;

/// Checks if the i-th bit is set
#[inline]
pub fn is_bit_set(x: &BigUint, i: u64) -> bool {
    get_bit(x, i) == 1
}

/// Returns the i-th bit.
#[inline]
fn get_bit(x: &BigUint, i: u64) -> u8 {
    let j = (i / BITS) as u64;
    // if is out of range of the set words, it is always false.
    let bits: u64 = x.bits();
    if i >= bits {
        return 0;
    }

    (x.get_limb(j) >> (i % BITS) & 1) as u8
}

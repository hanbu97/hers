use rug::Integer;

pub mod barrett;
pub mod conditional;
pub mod montgomery;

/// Helper function to compute high and low 64 bits of 128-bit product
#[inline(always)]
pub fn mul_hi_lo(a: u64, b: u64) -> (u64, u64) {
    let result = (a as u128) * (b as u128);
    ((result >> 64) as u64, result as u64)
}

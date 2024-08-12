// ref: https://github.com/dignifiedquire/num-bigint
use num_bigint::{BigInt, BigUint};
use num_bigint_dig::algorithms::jacobi;
use num_bigint_dig::RandBigInt;
use num_integer::Integer;
use num_traits::FromPrimitive;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

pub trait RandBigInt {
    /// Generate a random `BigUint` of the given bit size.
    fn gen_biguint(&mut self, bit_size: u64) -> BigUint;

    /// Generate a random `BigUint` less than the given bound. Fails
    /// when the bound is zero.
    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint;

    /// Generate a random `BigUint` within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint;

    /// Generate a random `BigInt` within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt;
}

impl<R: Rng + ?Sized> RandBigInt for R {
    fn gen_biguint(&mut self, bit_size: u64) -> BigUint {
        if bit_size == 0 {
            return BigUint::zero();
        }

        let num_u32 = ((bit_size + 31) / 32) as usize;
        let mut data = Vec::with_capacity(num_u32);

        // Generate random u32 values
        for _ in 0..num_u32 - 1 {
            data.push(self.gen::<u32>());
        }

        // Handle the last (possibly partial) u32
        let rem_bits = bit_size % 32;
        let last_u32 = if rem_bits == 0 {
            self.gen::<u32>()
        } else {
            self.gen::<u32>() & ((1 << rem_bits) - 1)
        };
        data.push(last_u32);

        // Create BigUint from the generated u32 slice
        let mut result = BigUint::from_slice(&data);

        // Ensure the most significant bit is set
        if !result.bit((bit_size - 1) as u64) {
            result.set_bit((bit_size - 1) as u64, true);
        }

        result
    }

    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint {
        assert!(!bound.is_zero());
        let bits = bound.bits();
        loop {
            let n = self.gen_biguint(bits);
            if n < *bound {
                return n;
            }
        }
    }

    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            self.gen_biguint_below(ubound)
        } else {
            lbound + self.gen_biguint_below(&(ubound - lbound))
        }
    }

    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            BigInt::from(self.gen_biguint_below(ubound.magnitude()))
        } else if ubound.is_zero() {
            lbound + BigInt::from(self.gen_biguint_below(lbound.magnitude()))
        } else {
            let delta = ubound - lbound;
            lbound + BigInt::from(self.gen_biguint_below(delta.magnitude()))
        }
    }
}

// Define a new trait for extending BigUint functionality
pub trait BigUintExt {
    fn get_limb(&self, i: u64) -> u32;
}

// Implement the new trait for BigUint
impl BigUintExt for BigUint {
    fn get_limb(&self, i: u64) -> u32 {
        // Convert BigUint to a slice of u32
        let digits = self.to_u32_digits();

        // Return 0 if the index is out of range, otherwise return the limb
        if i >= digits.len() as u64 {
            0
        } else {
            digits[i as usize]
        }
    }
}

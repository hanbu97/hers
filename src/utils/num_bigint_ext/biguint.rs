use super::*;

/// A generic trait for converting a value to a `BigUint`.
pub trait ToBigUint {
    /// Converts the value of `self` to a `BigUint`.
    fn to_biguint(&self) -> Option<BigUint>;
}

impl ToBigUint for BigUint {
    #[inline]
    fn to_biguint(&self) -> Option<BigUint> {
        Some(self.clone())
    }
}

/// A generic trait for converting a value to a `BigUint`, and consuming the value.
pub trait IntoBigUint {
    /// Converts the value of `self` to a `BigUint`.
    fn into_biguint(self) -> Option<BigUint>;
}

impl IntoBigUint for BigUint {
    #[inline]
    fn into_biguint(self) -> Option<BigUint> {
        Some(self)
    }
}

macro_rules! impl_to_biguint {
    ($T:ty, $from_ty:path) => {
        impl ToBigUint for $T {
            #[inline]
            fn to_biguint(&self) -> Option<BigUint> {
                $from_ty(*self)
            }
        }

        impl IntoBigUint for $T {
            #[inline]
            fn into_biguint(self) -> Option<BigUint> {
                $from_ty(self)
            }
        }
    };
}

impl_to_biguint!(isize, FromPrimitive::from_isize);
impl_to_biguint!(i8, FromPrimitive::from_i8);
impl_to_biguint!(i16, FromPrimitive::from_i16);
impl_to_biguint!(i32, FromPrimitive::from_i32);
impl_to_biguint!(i64, FromPrimitive::from_i64);
#[cfg(has_i128)]
impl_to_biguint!(i128, FromPrimitive::from_i128);

impl_to_biguint!(usize, FromPrimitive::from_usize);
impl_to_biguint!(u8, FromPrimitive::from_u8);
impl_to_biguint!(u16, FromPrimitive::from_u16);
impl_to_biguint!(u32, FromPrimitive::from_u32);
impl_to_biguint!(u64, FromPrimitive::from_u64);
#[cfg(has_i128)]
impl_to_biguint!(u128, FromPrimitive::from_u128);

impl_to_biguint!(f32, FromPrimitive::from_f32);
impl_to_biguint!(f64, FromPrimitive::from_f64);

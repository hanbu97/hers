use super::*;

use num_traits::CheckedSub;

#[must_use]
#[inline]
pub(crate) fn distance(x: &BigUint, y: &BigUint) -> BigUint {
    x.checked_sub(y).unwrap_or_else(|| y - x)
}

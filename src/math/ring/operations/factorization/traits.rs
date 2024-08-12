use super::*;

pub trait CompositeSplitter {
    /// Undefined behavior if `n` is prime.
    fn divisor(&self, n: &BigUint) -> BigUint;

    fn split(&self, n: &BigUint) -> (BigUint, BigUint) {
        let d1 = self.divisor(n);
        let d2 = n / &d1;
        if d1 < d2 {
            (d1, d2)
        } else {
            (d2, d1)
        }
    }
}

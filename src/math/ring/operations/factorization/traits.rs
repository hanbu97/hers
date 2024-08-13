use super::*;

pub trait CompositeSplitter {
    /// Undefined behavior if `n` is prime.
    fn divisor(&self, n: &Integer) -> Integer;

    fn split(&self, n: &Integer) -> (Integer, Integer) {
        let d1 = self.divisor(n);
        let d2: Integer = (n / &d1).into();
        if d1 < d2 {
            (d1, d2)
        } else {
            (d2, d1)
        }
    }
}

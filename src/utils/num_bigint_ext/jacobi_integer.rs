use num_integer::Integer as _;
use num_traits::identities::One;
use rug::Integer;

// Define a new trait for extending BigUint functionality
pub trait IntegerExt {
    fn get_limb(&self, i: u64) -> u32;
}

impl IntegerExt for rug::Integer {
    fn get_limb(&self, i: u64) -> u32 {
        let limbs = self.as_limbs();
        if i >= limbs.len() as u64 {
            0
        } else {
            limbs[i as usize] as u32
        }
    }
}

/// Jacobi returns the Jacobi symbol (x/y), either +1, -1, or 0.
/// The y argument must be an odd integer.
pub fn jacobi(x: &Integer, y: &Integer) -> isize {
    if !y.is_odd() {
        panic!(
            "invalid arguments, y must be an odd integer,but got {:?}",
            y
        );
    }

    let mut a = x.clone();
    let mut b = y.clone();
    let mut j = 1;

    if b.is_negative() {
        if a.is_negative() {
            j = -1;
        }
        b = -b;
    }

    loop {
        if b.is_one() {
            return j;
        }
        if a.is_zero() {
            return 0;
        }

        a = a.mod_floor(&b);
        if a.is_zero() {
            return 0;
        }

        // handle factors of 2 in a
        let s = a.find_one(0).unwrap();
        if s & 1 != 0 {
            let bmod8 = b.get_limb(0) & 7;
            if bmod8 == 3 || bmod8 == 5 {
                j = -j;
            }
        }

        let c: Integer = (&a >> s).into(); // a = 2^s*c

        // swap numerator and denominator
        if b.get_limb(0) & 3 == 3 && c.get_limb(0) & 3 == 3 {
            j = -j
        }

        a = b;
        b = c;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::FromPrimitive;

    #[test]
    fn test_jacobi() {
        let cases = [
            [0, 1, 1],
            [0, -1, 1],
            [1, 1, 1],
            [1, -1, 1],
            [0, 5, 0],
            [1, 5, 1],
            [2, 5, -1],
            [-2, 5, -1],
            [2, -5, -1],
            [-2, -5, 1],
            [3, 5, -1],
            [5, 5, 0],
            [-5, 5, 0],
            [6, 5, 1],
            [6, -5, 1],
            [-6, 5, 1],
            [-6, -5, -1],
        ];

        for case in cases.iter() {
            let x = Integer::from_i64(case[0]).unwrap();
            let y = Integer::from_i64(case[1]).unwrap();

            assert_eq!(case[2] as isize, jacobi(&x, &y), "jacobi({}, {})", x, y);
        }
    }
}

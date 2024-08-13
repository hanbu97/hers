use rug::{Assign, Integer};

use super::*;

// Implements Pollard's Rho algorithm for factorization using rug::Integer.
/// This function attempts to find a single factor of the input number.
pub fn get_factor_pollard_rho(m: &Integer) -> Integer {
    if m.is_prime() {
        return m.clone();
    }

    // Try different values of c in the polynomial x^2 + c
    for i in 1..10 {
        let mut x = Integer::from(2u32);
        let mut y = Integer::from(2u32);
        let mut d = Integer::from(1u32);
        let c = Integer::from(i);

        let one = Integer::from(1u32);

        while d.is_zero() || &d == m {
            // "Tortoise and hare" step
            x = polynomial_pollards_rho(&x, &c, m);
            y = polynomial_pollards_rho(&polynomial_pollards_rho(&y, &c, m), &c, m);

            d.assign(&x - &y);
            d = d.abs().gcd(m);

            if d != one {
                return d;
            }
        }
    }

    Integer::from(1u32) // Return 1 if no factor is found
}

/// Polynomial function used in Pollard's Rho algorithm: f(x) = x^2 + c mod n
fn polynomial_pollards_rho(x: &Integer, c: &Integer, n: &Integer) -> Integer {
    (x.clone().square() + c) % n
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_factor_pollard_rho() {
        let prime = 0x1fffffffffe00001u64;
        let m: Integer = Integer::from(prime - 1);
        assert_eq!(m.clone() % get_factor_pollard_rho(&m), Integer::from(0u32));
    }
}

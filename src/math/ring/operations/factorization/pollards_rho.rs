use super::*;

/// Implements Pollard's Rho algorithm for factorization.
/// This function attempts to find a single factor of the input number.
pub fn get_factor_pollard_rho(m: &BigUint) -> BigUint {
    if m.is_prime() {
        return m.clone();
    }

    // Try different values of c in the polynomial x^2 + c
    for i in 1..10 {
        let mut x = BigUint::from(2u32);
        let mut y = BigUint::from(2u32);
        let mut d = BigUint::one();
        let c = BigUint::from(i as u32);

        while d.is_zero() || &d == m {
            // "Tortoise and hare" step
            x = polynomial_pollards_rho(&x, &c, m);
            y = polynomial_pollards_rho(&polynomial_pollards_rho(&y, &c, m), &c, m);
            d = x.clone().sub(&y).gcd(m);
            if !d.is_one() {
                return d;
            }
        }
    }

    BigUint::one() // Return 1 if no factor is found
}

/// Polynomial function used in Pollard's Rho algorithm: f(x) = x^2 + c mod n
fn polynomial_pollards_rho(x: &BigUint, c: &BigUint, n: &BigUint) -> BigUint {
    (x.modpow(&BigUint::from(2u32), n) + c) % n
}

// ref: https://github.com/skyf0l/ecm-rs/blob/main/src/ecm.rs

pub mod ecm;
pub mod point;
// pub mod point_int;

use num_bigint::{BigUint, ToBigUint};
use num_integer::Integer as _;
use num_traits::One;

use crate::utils::{num_bigint_ext::rand::RandBigInt, prime::PrimeChecking};

use self::point::{CurveParams, Point};

/// Elliptic Curve Method (ECM) for factorization.
/// This function attempts to find a single factor of the input number.
pub fn get_factor_ecm(n: &BigUint) -> BigUint {
    let mut rng = rand::thread_rng();
    let two = 2.to_biguint().unwrap();
    let three = 3.to_biguint().unwrap();
    let four = 4.to_biguint().unwrap();
    let five = 5.to_biguint().unwrap();

    // Stage 1 parameters
    let b1: u32 = 10000; // Adjust as needed
    let b2: u32 = 100000; // Adjust as needed

    for _ in 0..100 {
        // Try up to 100 curves
        // Suyama's parameterization
        let sigma: BigUint = rng.gen_biguint_below(n);
        let u = (&sigma * &sigma - &five) % n;
        let v = (&four * &sigma) % n;

        let x = u.modpow(&three, n);
        let z = v.modpow(&three, n);

        let w = (&v - &u) % n;
        let a = ((&w * &w * &w * (&u + &u + &u + &v)) % n
            * Point::mod_inverse(&(&four * &x * &v), n).unwrap()
            - &two)
            % n;

        let a24 = (&a + &two) * Point::mod_inverse(&four, n).unwrap() % n;

        let params = CurveParams {
            a_24: a24,
            modulus: n.clone(),
        };

        let mut point = Point::new(x, z);

        // Stage 1
        for p in 2..=b1 {
            if p.is_prime() {
                let mut q = p;
                while q <= b1 {
                    point = point.mont_ladder(&q.to_biguint().unwrap(), &params);
                    q *= p;
                }
            }
        }

        let mut g = point.z_cord.gcd(n);
        if g > BigUint::one() && &g < n {
            return g;
        }

        // Stage 2 (simplified)
        for q in (b1 + 1)..=b2 {
            if q.is_prime() {
                point = point.mont_ladder(&q.to_biguint().unwrap(), &params);
                g = point.z_cord.gcd(n);
                if g > BigUint::one() && &g < n {
                    return g;
                }
            }
        }
    }

    n.clone() // Factor not found
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::{BigUint, ToBigUint};
    use num_traits::ToPrimitive;
    use num_traits::Zero;

    fn check_factor(n: &BigUint, factor: &BigUint) {
        assert!(n % factor == BigUint::zero());
        assert_ne!(factor, &BigUint::one());
        assert_ne!(factor, n);
    }

    #[test]
    fn test_small_composite() {
        let n = BigUint::from(1829u32);
        let factor = get_factor_ecm(&n);
        check_factor(&n, &factor);
    }

    #[test]
    fn test_medium_composite() {
        let n = BigUint::from(398883434337287u64);
        let factor = get_factor_ecm(&n);
        check_factor(&n, &factor);
    }

    #[test]
    fn test_large_composite() {
        let n: BigUint = "168541512131094651323".parse().unwrap();
        let factor = get_factor_ecm(&n);

        check_factor(&n, &factor);
    }

    #[test]
    fn test_very_large_composite() {
        let n: BigUint = "4269021180054189416198169786894227".parse().unwrap();
        let factor = get_factor_ecm(&n);

        check_factor(&n, &factor);
    }

    #[test]
    fn test_small_prime() {
        let n = BigUint::from(17u32);
        let factor = get_factor_ecm(&n);
        assert_eq!(factor, n);
    }

    #[test]
    fn test_large_prime() {
        let n: BigUint = "21472883178031195225853317139".parse().unwrap();
        let factor = get_factor_ecm(&n);
        assert_eq!(factor, n);
    }

    #[test]
    fn test_power_of_prime() {
        let n = BigUint::from(2u32).pow(16); // 65536
        let factor = get_factor_ecm(&n);
        assert_eq!(factor, BigUint::from(2u32));
    }

    #[test]
    fn test_product_of_primes() {
        let n =
            BigUint::from(2u32) * BigUint::from(3u32) * BigUint::from(5u32) * BigUint::from(7u32);
        let factor = get_factor_ecm(&n);
        check_factor(&n, &factor);
        assert!(vec![2u32, 3u32, 5u32, 7u32].contains(&factor.to_u32().unwrap()));
    }

    #[test]
    fn test_same_factors() {
        let n: BigUint = "7853316850129".parse().unwrap(); // 2802377 * 2802377
        let factor = get_factor_ecm(&n);
        check_factor(&n, &factor);
        assert_eq!(factor, BigUint::from(2802377u32));
    }
}

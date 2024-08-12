pub mod ecm;
pub mod pollards_rho;
pub mod traits;
pub mod utils;

use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_integer::Integer;

use num_traits::{One, Zero};
use rand::Rng;
use std::ops::Sub;

use crate::utils::{num_bigint_ext::rand::RandBigInt, prime::PrimeChecking};

use self::pollards_rho::get_factor_pollard_rho;

/// Finds all prime factors of a given BigUint.
/// Returns a sorted vector of prime factors.
pub fn get_factors(m: &BigUint) -> Vec<BigUint> {
    let mut m_cpy = m.clone();
    if m_cpy.is_prime() {
        return vec![m_cpy];
    }

    let mut f = std::collections::HashSet::new();

    // First, loop through small prime factors
    for &small_prime in SMALL_PRIMES.iter() {
        let small_prime = BigUint::from(small_prime);
        let mut add_factor = false;
        while (&m_cpy % &small_prime).is_zero() {
            m_cpy /= &small_prime;
            add_factor = true;
        }
        if add_factor {
            f.insert(small_prime);
        }
    }

    // Second, find the remaining large prime factors
    while !m_cpy.is_one() {
        if m_cpy.is_prime() {
            f.insert(m_cpy.clone());
            break;
        }

        // Try Pollard's Rho algorithm first
        let mut factor = get_factor_pollard_rho(&m_cpy);
        if factor.is_one() || factor == m_cpy {
            // If Pollard's Rho fails, try ECM factorization
            // factor = get_factor_ecm(&m_cpy);
        }

        // Remove all instances of this factor from m_cpy
        while (&m_cpy % &factor).is_zero() {
            m_cpy /= &factor;
        }

        f.insert(factor);
    }

    // Convert the set of factors to a sorted vector
    let mut factors: Vec<BigUint> = f.into_iter().collect();
    factors.sort();
    factors
}

/// Elliptic Curve Method (ECM) for factorization.
/// This function attempts to find a single factor of the input number.
// fn get_factor_ecm(n: &BigUint) -> BigUint {
//     if n.is_prime() {
//         return n.clone();
//     }

//     let mut ecm = ECM::new(n);
//     let one = BigUint::one();

//     loop {
//         // Generate a random elliptic curve and point
//         let (weierstrass, mut p) = Weierstrass::new_random_curve(n);
//         ecm.weierstrass = weierstrass;

//         let mut bound = 0.0;
//         let mut i = 2u64;

//         // Main ECM loop
//         while bound < ecm.b {
//             let (new_p, factor) = ecm.check_then_mul(i, &p);
//             if factor != one {
//                 return factor;
//             }
//             p = new_p;
//             i += 1;
//             bound += 1.0;
//         }
//     }
// }

// /// Struct representing the Elliptic Curve Method context
// struct ECM {
//     weierstrass: Weierstrass,
//     n: BigUint,
//     b: f64,
// }

// impl ECM {
//     /// Creates a new ECM context for the given number to factorize
//     fn new(n: &BigUint) -> Self {
//         let log_n = (n.bits() as f64 + 1.0) / 1.4426950408889634;
//         ECM {
//             weierstrass: Weierstrass::default(),
//             n: n.clone(),
//             b: ((2.0 * log_n * log_n.ln()).sqrt()).exp(),
//         }
//     }

//     /// Checks for a factor while adding two points on the elliptic curve
//     fn check_then_add(&self, p: &Point, q: &Point) -> (Point, BigUint) {
//         let n = &self.n;
//         let one = BigUint::one();

//         if p.x == q.x && p.y == q.y {
//             let gcd = (&p.y + &p.y).gcd(n);
//             if gcd != one {
//                 return (Point::default(), gcd);
//             }
//         } else {
//             let gcd = (&q.x - &p.x).gcd(n);
//             if gcd != one {
//                 return (Point::default(), gcd);
//             }
//         }

//         (self.weierstrass.add(p, q), one)
//     }

//     /// Checks for a factor while multiplying a point on the elliptic curve
//     fn check_then_mul(&self, k: u64, p: &Point) -> (Point, BigUint) {
//         let mut q = Point::zero();
//         let one = BigUint::one();

//         let mut k = k;
//         let mut p = p.clone();

//         // Double-and-add algorithm for point multiplication
//         while k > 0 {
//             if k & 1 == 1 {
//                 let (new_q, gcd) = self.check_then_add(&p, &q);
//                 if gcd != one {
//                     return (Point::default(), gcd);
//                 }
//                 q = new_q;
//             }

//             k >>= 1;

//             if k > 0 {
//                 let (new_p, gcd) = self.check_then_add(&p, &p);
//                 if gcd != one {
//                     return (Point::default(), gcd);
//                 }
//                 p = new_p;
//             }
//         }

//         (q, one)
//     }
// }

/// Represents a point on an elliptic curve
#[derive(Clone, Default)]
struct Point {
    x: BigUint,
    y: BigUint,
}

impl Point {
    /// Creates a point at infinity (identity element for the elliptic curve group)
    fn zero() -> Self {
        Point {
            x: BigUint::zero(),
            y: BigUint::one(),
        }
    }
}

// /// Represents a Weierstrass form elliptic curve: y^2 = x^3 + ax + b (mod p)
// #[derive(Default)]
// pub struct Weierstrass {
//     a: BigUint,
//     b: BigUint,
//     p: BigUint,
// }

// impl Weierstrass {
//     /// Generates a random elliptic curve and a point on it
//     fn new_random_curve(p: &BigUint) -> (Self, Point) {
//         let mut rng = rand::thread_rng();
//         // let x = rng.gen_range(BigUint::zero()..p.clone());
//         // let y = rng.gen_range(BigUint::zero()..p.clone());
//         // let a = rng.gen_range(BigUint::zero()..p.clone());

//         let x = rng.gen_biguint_range(&BigUint::zero(), p);
//         let y = rng.gen_biguint_range(&BigUint::zero(), p);
//         let a = rng.gen_biguint_range(&BigUint::zero(), p);

//         let x_cube = x.modpow(&BigUint::from(3u32), p);
//         let ax = (&a * &x) % p;
//         let y_square = (&y * &y) % p;
//         let b = ((&y_square - &x_cube) + &ax) % p;

//         (Weierstrass { a, b, p: p.clone() }, Point { x, y })
//     }

//     /// Adds two points on the elliptic curve
//     fn add(&self, p: &Point, q: &Point) -> Point {
//         if p.x == q.x && p.y == q.y {
//             return self.double(p);
//         }

//         let p = &self.p;
//         let m = ((&q.y - &p.y) * mod_inverse(&(&q.x - &p.x), p).unwrap()) % p;
//         let x = (&m * &m - &p.x - &q.x) % p;
//         let y = (&m * (&p.x - &x) - &p.y) % p;

//         Point { x, y }
//     }

//     /// Doubles a point on the elliptic curve
//     fn double(&self, p: &Point) -> Point {
//         let p_mod = &self.p;
//         let m = ((BigUint::from(3u32) * &p.x * &p.x + &self.a)
//             * mod_inverse(&(BigUint::from(2u32) * &p.y), p_mod).unwrap())
//             % p_mod;
//         let x = (&m * &m - BigUint::from(2u32) * &p.x) % p_mod;
//         let y = (&m * (&p.x - &x) - &p.y) % p_mod;

//         Point { x, y }
//     }
// }

// /// Computes the modular inverse of a modulo m
// fn mod_inverse(a: &BigUint, m: &BigUint) -> Option<BigUint> {
//     let (g, x, _) = extended_gcd(a, m);
//     if g != BigUint::one() {
//         None
//     } else {
//         Some((x % m + m) % m)
//     }
// }

// /// Computes the extended Euclidean algorithm
// fn extended_gcd(a: &BigUint, b: &BigUint) -> (BigUint, BigInt, BigInt) {
//     if b.is_zero() {
//         (a.clone(), BigInt::one(), BigInt::zero())
//     } else {
//         let (g, x, y) = extended_gcd(b, &(a % b));
//         (g, y, x - (a / b * y))
//     }
// }

/// A list of small prime numbers used for trial division
const SMALL_PRIMES: &[u64] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523,
    541,
    // ... (其余的素数被省略以节省空间)
];

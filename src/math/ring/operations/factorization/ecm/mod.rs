// ref: https://github.com/skyf0l/ecm-rs/blob/main/src/ecm.rs

pub mod point;
// pub mod point_int;

use num_bigint::{BigUint, ToBigUint};
use num_integer::Integer as _;
use num_traits::{One, Zero};

use crate::utils::{num_bigint_ext::rand::RandBigInt, prime::PrimeChecking};

#[derive(Clone, Debug)]
struct Point {
    x: BigUint,
    z: BigUint,
}

impl Point {
    fn new(x: BigUint, z: BigUint) -> Self {
        Point { x, z }
    }
}

/// Modular exponentiation: (base^exponent) % modulus
fn mod_pow(base: &BigUint, exponent: &BigUint, modulus: &BigUint) -> BigUint {
    let mut result = BigUint::one();
    let mut base = base.clone();
    let mut exp = exponent.clone();

    while !exp.is_zero() {
        if exp.is_odd() {
            result = (result * &base) % modulus;
        }
        exp >>= 1;
        base = (&base * &base) % modulus;
    }
    result
}

/// Modular inverse: returns a such that (a * n) % modulus == 1
fn mod_inverse(n: &BigUint, modulus: &BigUint) -> Option<BigUint> {
    let (mut t, mut newt) = (BigUint::zero(), BigUint::one());
    let (mut r, mut newr) = (modulus.clone(), n.clone());

    while !newr.is_zero() {
        let quotient = &r / &newr;
        (t, newt) = (newt.clone(), t - quotient.clone() * newt);
        (r, newr) = (newr.clone(), r - quotient * newr);
    }

    if r > BigUint::one() {
        None
    } else if t < BigUint::zero() {
        Some(t + modulus)
    } else {
        Some(t)
    }
}

/// Addition of points on the elliptic curve
fn add_points(p: &Point, q: &Point, diff: &Point, a24: &BigUint, n: &BigUint) -> Point {
    let u = (&p.x - &p.z) * (&q.x + &q.z) % n;
    let v = (&p.x + &p.z) * (&q.x - &q.z) % n;
    let x = (&diff.z * (&u + &v).pow(2u32)) % n;
    let z = (&diff.x * (&u - &v).pow(2u32)) % n;
    Point::new(x, z)
}

/// Doubling a point on the elliptic curve
fn double_point(p: &Point, a24: &BigUint, n: &BigUint) -> Point {
    let u = (&p.x + &p.z).pow(2u32) % n;
    let v = (&p.x - &p.z).pow(2u32) % n;
    let x = (&u * &v) % n;
    let t = u - &v;
    let z = (t.clone() * (a24 * &t + &v)) % n;
    Point::new(x, z)
}

/// Montgomery ladder for scalar multiplication
fn montgomery_ladder(k: &BigUint, base: &Point, a24: &BigUint, n: &BigUint) -> Point {
    let mut r0 = base.clone();
    let mut r1 = double_point(base, a24, n);

    for i in (0..k.bits()).rev() {
        if k.bit(i) {
            r0 = add_points(&r0, &r1, base, a24, n);
            r1 = double_point(&r1, a24, n);
        } else {
            r1 = add_points(&r0, &r1, base, a24, n);
            r0 = double_point(&r0, a24, n);
        }
    }
    r0
}

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

        let x = mod_pow(&u, &three, n) % n;
        let z = mod_pow(&v, &three, n) % n;

        let w = (&v - &u) % n;
        let a = ((&w * &w * &w * (&u + &u + &u + &v)) % n
            * mod_inverse(&(&four * &x * &v), n).unwrap()
            - &two)
            % n;

        let a24 = (&a + &two) * mod_inverse(&four, n).unwrap() % n;

        let mut point = Point::new(x, z);

        // Stage 1
        for p in 2..=b1 {
            if p.is_prime() {
                let mut q = p;
                while q <= b1 {
                    point = montgomery_ladder(&q.to_biguint().unwrap(), &point, &a24, n);
                    q *= p;
                }
            }
        }

        let mut g = point.z.gcd(n);
        if g > BigUint::one() && &g < n {
            return g;
        }

        // Stage 2 (simplified)
        for q in (b1 + 1)..=b2 {
            if q.is_prime() {
                point = montgomery_ladder(&q.to_biguint().unwrap(), &point, &a24, n);
                g = point.z.gcd(n);
                if g > BigUint::one() && &g < n {
                    return g;
                }
            }
        }
    }

    n.clone() // Factor not found
}

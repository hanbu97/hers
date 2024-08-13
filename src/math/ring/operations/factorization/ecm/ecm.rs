use crate::utils::num_bigint_ext::rand::RandBigInt;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use primal::Primes;
use std::collections::HashMap;

use crate::utils::prime::PrimeChecking;
use num_integer::Integer;
use rand::SeedableRng;

use super::point::{CurveParams, Point};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("Bounds should be an even integer")]
    BoundsNotEven,
    #[error("Too small bounds")]
    BoundsTooSmall,
    #[error("The factorization failed")]
    ECMFailed,
    #[error("The number is prime")]
    NumberIsPrime,
}

pub fn ecm_one_factor(
    n: &BigUint,
    b1: usize,
    b2: usize,
    max_curve: usize,
    rng: &mut impl rand::Rng,
) -> Result<BigUint, Error> {
    if b1 % 2 != 0 || b2 % 2 != 0 {
        return Err(Error::BoundsNotEven);
    }

    if n.is_prime() {
        return Err(Error::NumberIsPrime);
    }

    let mut curve = 0;
    let d = (b2 as f64).sqrt() as usize;
    let two_d = 2 * d;
    let mut beta: Vec<BigUint> = vec![BigUint::zero(); d + 1];
    let mut s: Vec<Point> = vec![Point::new(BigUint::zero(), BigUint::zero()); d + 1];
    let mut k = BigUint::one();

    for p in Primes::all().take_while(|&p| p <= b1) {
        let p = BigUint::from(p);
        let mut power = p.clone();
        while power <= BigUint::from(b1) {
            k *= &p;
            power *= &p;
        }
    }

    while curve <= max_curve {
        curve += 1;

        // Suyama's Parametrization
        let sigma: BigUint = rng.gen_biguint_below(n);

        // dbg!(&(&sigma * &sigma - BigUint::from(5u32)) % n);
        // let u = (&sigma * &sigma - BigUint::from(5u32)) % n;
        let u = (&sigma * &sigma + n - BigUint::from(5u32)) % n;
        let v = (BigUint::from(4u32) * &sigma) % n;

        let diff = (&v + n - &u) % n;

        let u_3 = u.modpow(&BigUint::from(3u32), n);
        let v_3 = v.modpow(&BigUint::from(3u32), n);

        let c = match (BigUint::from(4u32) * &u_3 * &v).modinv(n) {
            Some(c) => {
                (diff.modpow(&BigUint::from(3u32), n) * (BigUint::from(4u32) * &u + &v) * c
                    - BigUint::from(2u32))
                    % n
            }
            None => return Ok((BigUint::from(4u32) * u_3 * v).gcd(n)),
        };

        let a24 = (&c + BigUint::from(2u32)) * BigUint::from(4u32).modinv(n).unwrap() % n;
        let params = CurveParams {
            a_24: a24,
            modulus: n.clone(),
        };
        let q = Point::new(u_3, v_3);
        let q = q.mont_ladder(&k, &params);
        let g = q.z_cord.gcd(n);

        // Stage 1 factor
        if &g != n && g != BigUint::one() {
            return Ok(g);
        }

        // Stage 1 failure. Q.z = 0, Try another curve
        if &g == n {
            continue;
        }

        // Stage 2 - Improved Standard Continuation
        s[1] = q.double(&params);
        s[2] = s[1].double(&params);
        beta[1] = (&s[1].x_cord * &s[1].z_cord) % n;
        beta[2] = (&s[2].x_cord * &s[2].z_cord) % n;

        for d in 3..=d {
            s[d] = s[d - 1].add(&s[1], &s[d - 2], &params);
            beta[d] = (&s[d].x_cord * &s[d].z_cord) % n;
        }

        let mut g = BigUint::one();
        let b = b1 - 1;
        let mut t = q.mont_ladder(&BigUint::from(b - two_d), &params);
        let mut r = q.mont_ladder(&BigUint::from(b), &params);

        let mut primes = Primes::all().skip_while(|&q| q < b);
        for rr in (b..b2).step_by(two_d) {
            let alpha = (&r.x_cord * &r.z_cord) % n;
            for q in primes.by_ref().take_while(|&q| q <= rr + two_d) {
                let delta = (q - rr) / 2;

                // dbg!("************");
                // dbg!((&r.x_cord - &s[d].x_cord));
                // dbg!("************");

                // let f =
                //     (&r.x_cord - &s[d].x_cord) * (&r.z_cord + &s[d].z_cord) - &alpha + &beta[delta];

                let f = ((&r.x_cord + n - &s[d].x_cord) % n * (&r.z_cord + &s[d].z_cord) % n + n
                    - &alpha
                    + &beta[delta])
                    % n;

                g = (g * f) % n;
            }
            // Swap
            std::mem::swap(&mut t, &mut r);
            r = r.add(&s[d], &t, &params);
        }
        g = g.gcd(n);

        // Stage 2 Factor found
        if &g != n && g != BigUint::one() {
            return Ok(g);
        }
    }

    // ECM failed, Increase the bounds
    Err(Error::ECMFailed)
}

fn optimal_params(digits: usize) -> (usize, usize, usize) {
    match digits {
        1..=10 => (2_000, 160_000, 35),
        11..=15 => (5_000, 500_000, 500),
        16..=20 => (11_000, 1_900_000, 74),
        21..=25 => (50_000, 13_000_000, 214),
        26..=30 => (250_000, 130_000_000, 430),
        31..=35 => (1_000_000, 1_000_000_000, 904),
        36..=40 => (3_000_000, 5_700_000_000, 2350),
        41..=45 => (11_000_000, 35_000_000_000, 4480),
        46..=50 => (44_000_000, 240_000_000_000, 7553),
        51..=55 => (110_000_000, 780_000_000_000, 17769),
        56..=60 => (260_000_000, 3_200_000_000_000, 42017),
        _ => (850_000_000, 16_000_000_000_000, 69408),
    }
}

pub fn ecm(n: &BigUint) -> Result<HashMap<BigUint, usize>, Error> {
    let optimal_params = optimal_params(n.to_str_radix(10).len());
    ecm_with_params(
        n,
        optimal_params.0,
        optimal_params.1,
        optimal_params.2,
        1234,
    )
}

pub fn ecm_with_params(
    n: &BigUint,
    b1: usize,
    b2: usize,
    max_curve: usize,
    seed: u64,
) -> Result<HashMap<BigUint, usize>, Error> {
    let mut factors = HashMap::new();
    let mut n = n.clone();

    for prime in Primes::all().take(100_000) {
        let prime = BigUint::from(prime);
        while (&n % &prime).is_zero() {
            n /= &prime;
            *factors.entry(prime.clone()).or_insert(0) += 1;
        }
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    while n != BigUint::one() {
        let factor = ecm_one_factor(&n, b1, b2, max_curve, &mut rng).unwrap_or_else(|_| n.clone());

        while (&n % &factor).is_zero() {
            n /= &factor;
            *factors.entry(factor.clone()).or_insert(0) += 1;
        }
    }

    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::{BigUint, ToBigUint};

    fn ecm(n: &BigUint) -> Result<HashMap<BigUint, usize>, Error> {
        super::ecm(n)
    }

    #[test]
    fn sympy_1() {
        let n = 398883434337287u64.to_biguint().unwrap();
        let factors = ecm(&n).unwrap();
        assert_eq!(factors.len(), 2);
        assert!(factors.contains_key(&99476569u64.to_biguint().unwrap()));
        assert!(factors.contains_key(&4009823u64.to_biguint().unwrap()));
    }

    #[test]
    fn sympy_2() {
        let n = 46167045131415113u64.to_biguint().unwrap();
        let factors = ecm(&n).unwrap();
        assert_eq!(factors.len(), 3);
        assert!(factors.contains_key(&43u64.to_biguint().unwrap()));
        assert!(factors.contains_key(&2634823u64.to_biguint().unwrap()));
        assert!(factors.contains_key(&407485517u64.to_biguint().unwrap()));
    }

    #[test]
    fn sympy_3() {
        let n = 64211816600515193u64.to_biguint().unwrap();
        let factors = ecm(&n).unwrap();
        assert_eq!(factors.len(), 3);
        assert!(factors.contains_key(&281719u64.to_biguint().unwrap()));
        assert!(factors.contains_key(&359641u64.to_biguint().unwrap()));
        assert!(factors.contains_key(&633767u64.to_biguint().unwrap()));
    }

    #[test]
    fn same_factors() {
        let n = 7853316850129u64.to_biguint().unwrap();
        let factors = ecm(&n).unwrap();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors[&2802377u64.to_biguint().unwrap()], 2);
    }

    #[test]
    fn small_prime() {
        let n = 17u64.to_biguint().unwrap();
        let result = ecm(&n);
        assert!(matches!(result, Err(Error::NumberIsPrime)));
    }

    #[test]
    fn big_prime() {
        let n = "21472883178031195225853317139".parse::<BigUint>().unwrap();
        let result = ecm(&n);
        assert!(matches!(result, Err(Error::NumberIsPrime)));
    }
}

use std::collections::HashMap;

use num_bigint::BigUint;

// fn get_factor_pollard_rho(m: &BigInt) -> Option<BigInt> {
//     if is_prime(m) {
//         return Some(m.clone());
//     }

//     for c in 1..10 {
//         let mut x = BigInt::from(2u32);
//         let mut y = BigInt::from(2u32);
//         let mut d = BigInt::one();

//         while d.is_one() {
//             x = polynomial_pollards_rho(&x, &c.to_bigint().unwrap(), m);
//             y = polynomial_pollards_rho(
//                 &polynomial_pollards_rho(&y, &c.to_bigint().unwrap(), m),
//                 &c.to_bigint().unwrap(),
//                 m,
//             );
//             d = (&x - &y).abs().gcd(m);

//             if d != BigInt::one() && d != m.clone() {
//                 return Some(d);
//             }
//         }
//     }

//     None
// }

// pub fn get_factors(m: &BigUint) -> Vec<BigUint> {
//     let mut m_copy = m.clone();
//     let mut factors = HashMap::new();

//     if m_copy.is_prime() {
//         return vec![m_copy];
//     }

//     // Check small prime factors
//     for &small_prime in SMALL_PRIMES.iter() {
//         let small_prime = small_prime.to_bigint().unwrap();
//         while (&m_copy % &small_prime).is_zero() {
//             m_copy /= &small_prime;
//             factors.entry(small_prime.clone()).or_insert(0) += 1;
//         }
//     }

//     // Find remaining large prime factors
//     while !m_copy.is_one() {
//         if is_prime(&m_copy) {
//             *factors.entry(m_copy.clone()).or_insert(0) += 1;
//             break;
//         }

//         let factor = if let Some(f) = get_factor_pollard_rho(&m_copy) {
//             f
//         } else {
//             get_factor_ecm(&m_copy)
//         };

//         while (&m_copy % &factor).is_zero() {
//             m_copy /= &factor;
//         }
//         *factors.entry(factor).or_insert(0) += 1;
//     }

//     factors.into_iter().map(|(k, _)| k).collect()
// }

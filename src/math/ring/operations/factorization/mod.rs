pub mod ecm;
pub mod pollards_rho;
pub mod traits;
use num_traits::One;

use rug::Integer;

use crate::utils::prime::PrimeChecking;

use self::{ecm::get_factor_ecm, pollards_rho::get_factor_pollard_rho};

/// Finds all prime factors of a given BigUint.
/// Returns a sorted vector of prime factors.
pub fn get_factors(m: &Integer) -> anyhow::Result<Vec<Integer>> {
    let mut m_cpy = m.clone();
    if m_cpy.is_prime() {
        return Ok(vec![m_cpy]);
    }

    let mut f = std::collections::HashSet::new();

    for prime in primal::Primes::all().take(10000) {
        let small_prime = Integer::from(prime);
        let mut add_factor = false;

        while m_cpy.is_divisible_u(prime as u32) {
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
        let mut factor: Integer = get_factor_pollard_rho(&m_cpy);
        if factor.is_one() || factor == m_cpy {
            // If Pollard's Rho fails, try ECM factorization
            factor = get_factor_ecm(&m_cpy)?;
        }

        // Remove all instances of this factor from m_cpy
        let temp: Integer = (&m_cpy % &factor).into();
        while temp.is_zero() {
            m_cpy /= &factor;
        }

        f.insert(factor);
    }

    // Convert the set of factors to a sorted vector
    let mut factors: Vec<Integer> = f.into_iter().collect();
    factors.sort();
    Ok(factors)
}

/// Checks if the given factors completely factorize the input number.
///
/// # Arguments
///
/// * `p` - The number to be factorized.
/// * `factors` - A slice of factors to check against.
///
/// # Returns
///
/// `true` if the factors completely factorize `p`, `false` otherwise.
pub fn check_factorization(p: &Integer, factors: &[Integer]) -> bool {
    let mut remaining = p.clone();
    let zero = Integer::from(0);

    for factor in factors {
        while (remaining.clone() % factor) == zero {
            remaining /= factor;
        }
    }

    remaining.is_one()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_small_primes() {
        let max = primal::Primes::all().take(10000).last().unwrap();
        dbg!(max);
    }

    #[test]
    fn test_get_factors() {
        let m = Integer::from(0x1fffffffffe00001u64 - 1);
        let factors = get_factors(&m).unwrap();

        println!("m: {} {:?}", m, factors);
        assert!(check_factorization(&m, &factors));
    }
}

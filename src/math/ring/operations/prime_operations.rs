use rug::Integer;

use super::{factorization::get_factors, mod_exp};
use crate::{math::ring::errors::MathError, utils::prime::PrimeChecking};

/// Computes the smallest primitive root of the given prime q.
/// The unique factors of q-1 can be given to speed up the search for the root.
pub fn primitive_root(q: u64, factors: Option<Vec<u64>>) -> anyhow::Result<(u64, Vec<u64>)> {
    let factors = match factors {
        Some(f) => {
            check_factors(q - 1, &f)?;
            f
        }
        None => {
            let factors_big = get_factors(&Integer::from(q - 1))?;
            factors_big.iter().map(|f| f.to_u64().unwrap()).collect()
        }
    };

    let mut g = 2u64; // g = 3 ?
    loop {
        let mut is_primitive_root = true;
        for &factor in &factors {
            if mod_exp(g, (q - 1) / factor, q) == 1 {
                is_primitive_root = false;
                break;
            }
        }
        if is_primitive_root {
            return Ok((g, factors));
        }
        g += 1;
    }
}

/// Checks that the given list of factors contains all the unique primes of m.
pub fn check_factors(m: u64, factors: &[u64]) -> Result<(), MathError> {
    let mut m = m;
    for &factor in factors {
        if !factor.is_prime() {
            return Err(MathError::CompositeFactor);
        }
        while m % factor == 0 {
            m /= factor;
        }
    }
    if m != 1 {
        Err(MathError::IncompleteFactorList)
    } else {
        Ok(())
    }
}

/// Checks that g is a valid primitive root mod q, given the factors of q-1.
pub fn check_primitive_root(g: u64, q: u64, factors: &[u64]) -> Result<(), MathError> {
    check_factors(q - 1, factors)?;

    for &factor in factors {
        if mod_exp(g, (q - 1) / factor, q) == 1 {
            return Err(MathError::InvalidPrimitiveRoot);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_factors() {
        let test_cases = [
            (12, vec![2, 3], Ok(())),
            (30, vec![2, 3, 5], Ok(())),
            (12, vec![2, 2, 3], Ok(())),
            (12, vec![2], Err(MathError::IncompleteFactorList)),
            (17, vec![17], Ok(())),
            (100, vec![2, 5], Ok(())),
            (100, vec![2, 2, 5, 5], Ok(())),
            (100, vec![2, 5, 7], Ok(())),
            (1, vec![], Ok(())),
            (0, vec![], Err(MathError::IncompleteFactorList)),
            (2, vec![2], Ok(())),
            (
                4611686018427387904,
                vec![2, 4611686018427387903],
                Err(MathError::CompositeFactor),
            ),
            (
                18446744073709551615,
                vec![3, 6148914691236517205],
                Err(MathError::CompositeFactor),
            ),
            (
                18446744073709551614,
                vec![2, 9223372036854775807],
                Err(MathError::CompositeFactor),
            ),
        ];

        for (i, &(m, ref factors, ref expected_result)) in test_cases.iter().enumerate() {
            let result = check_factors(m, factors);
            // println!("Test case {}:", i + 1);
            // println!("  Input: m = {}, factors = {:?}", m, factors);
            // println!("  Output: {:?}\n", result);
            assert_eq!(result, *expected_result, "Test case {} failed", i + 1);
        }
    }

    #[test]
    fn test_check_primitive_root() {
        let test_cases = [
            (3, 7, vec![2, 3], Ok(())),
            (2, 7, vec![2, 3], Err(MathError::InvalidPrimitiveRoot)),
            (2, 11, vec![2, 5], Ok(())),
            (3, 11, vec![2, 5], Err(MathError::InvalidPrimitiveRoot)),
            (5, 23, vec![2, 11], Ok(())),
            (2, 23, vec![2, 11], Err(MathError::InvalidPrimitiveRoot)),
            (2, 37, vec![2, 3], Ok(())),
            (5, 61, vec![2, 3, 5], Err(MathError::InvalidPrimitiveRoot)),
            (2, 61, vec![2, 3, 5], Ok(())),
            (2, 4, vec![2], Err(MathError::IncompleteFactorList)),
            (3, 65537, vec![2], Ok(())),
            (5, 65537, vec![2], Ok(())),
            (2, 65537, vec![2], Err(MathError::InvalidPrimitiveRoot)),
        ];

        for (i, &(g, q, ref factors, ref expected_result)) in test_cases.iter().enumerate() {
            let result = check_primitive_root(g, q, factors);
            // println!("Test case {}:", i + 1);
            // println!("  Input: g = {}, q = {}, factors = {:?}", g, q, factors);
            // println!("  Output: {:?}\n", result);
            assert_eq!(result, *expected_result, "Test case {} failed", i + 1);
        }
    }
}

// ref: https://github.com/skyf0l/ecm-rs/blob/main/src/ecm.rs

use rand::{rngs::OsRng, Rng};
use rug::{rand::RandState, Integer};

use self::ecm::ecm_one_factor;

pub mod ecm;
pub mod errors;
pub mod point;

/// Elliptic Curve Method (ECM) for factorization.
/// This function attempts to find a single factor of the input number.
pub fn get_factor_ecm(n: &Integer) -> anyhow::Result<Integer> {
    let digits = n.to_string().len();
    let (b1, b2, max_curve) = ecm::optimal_params(digits);

    let mut rgen = RandState::new();
    // rgen.seed(&Integer::from(1234)); // You can change this seed if needed
    let mut csprng = OsRng;
    rgen.seed(&Integer::from(csprng.gen::<u64>()));

    match ecm_one_factor(n, b1, b2, max_curve, &mut rgen) {
        Ok(factor) => Ok(factor),
        Err(errors::ECMErrors::NumberIsPrime) => Ok(n.clone()),
        Err(e) => Err(e.into()),
    }
}

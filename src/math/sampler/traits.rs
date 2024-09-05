use std::fmt::Debug;

use crate::math::ring::polynomial::Poly;

/// Trait for polynomial samplers.
pub trait Sampler {
    /// Samples a polynomial and stores it in the given polynomial.
    fn read(&mut self, pol: &mut Poly);

    /// Samples a new polynomial and returns it.
    fn read_new(&mut self) -> Poly;

    /// Samples a polynomial and adds it to the given polynomial.
    fn read_and_add(&mut self, pol: &mut Poly);

    /// Returns a new sampler instance at the specified level.
    fn at_level(&self, level: usize) -> Self;
}

/// Trait for distribution parameters.
pub trait DistributionParameters: Debug {
    /// Returns the type of the distribution as a string.
    fn type_name(&self) -> &'static str;
}

/// Parameters for discrete Gaussian distribution.
#[derive(Debug, Clone, Copy)]
pub struct DiscreteGaussian {
    /// Standard deviation of the distribution.
    pub sigma: f64,
    /// Upper bound for the distribution.
    pub bound: f64,
}

impl DistributionParameters for DiscreteGaussian {
    fn type_name(&self) -> &'static str {
        "DiscreteGaussian"
    }
}

/// Parameters for ternary distribution.
#[derive(Debug, Clone, Copy)]
pub struct Ternary {
    /// Probability parameter for ternary distribution.
    pub p: f64,
    /// Hamming weight for sparse ternary distribution.
    pub h: usize,
}

impl DistributionParameters for Ternary {
    fn type_name(&self) -> &'static str {
        "Ternary"
    }
}

/// Parameters for uniform distribution.
#[derive(Debug, Clone, Copy)]
pub struct Uniform;

impl DistributionParameters for Uniform {
    fn type_name(&self) -> &'static str {
        "Uniform"
    }
}

/// Trait for pseudo-random number generators.
pub trait PRNG: std::io::Read {
    /// Resets the PRNG to its initial state.
    fn reset(&mut self);
}

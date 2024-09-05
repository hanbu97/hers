use std::fmt::Debug;

/// Represents a ring structure
pub trait Ring: Clone + Debug {
    /// The polynomial type associated with this ring
    type Poly: Clone + Debug;

    /// Creates a new ring with the given parameters
    fn new(n: usize, moduli: Vec<u64>) -> Self;

    /// Returns the degree of the ring
    fn n(&self) -> usize;

    /// Returns the modulus for the given level
    fn modulus(&self, level: usize) -> u64;

    /// Returns the current level of the ring
    fn level(&self) -> usize;

    /// Returns the mask for the given level
    fn mask(&self, level: usize) -> u64;

    /// Creates a new polynomial in this ring
    fn new_poly(&self) -> Self::Poly;

    /// Returns a new ring instance at the specified level
    fn at_level(&self, level: usize) -> Self;
}

/// Represents a polynomial sampler
pub trait Sampler {
    /// The ring type associated with this sampler
    type Ring: Ring;

    /// Samples a polynomial and stores it in the given polynomial
    fn read(&mut self, pol: &mut <Self::Ring as Ring>::Poly);

    /// Samples a new polynomial and returns it
    fn read_new(&mut self) -> <Self::Ring as Ring>::Poly;

    /// Samples a polynomial and adds it to the given polynomial
    fn read_and_add(&mut self, pol: &mut <Self::Ring as Ring>::Poly);

    /// Returns a new sampler instance at the specified level
    fn at_level(&self, level: usize) -> Self;
}

/// Represents distribution parameters
pub trait DistributionParameters: Debug {
    /// Returns the type of the distribution as a string
    fn type_name(&self) -> &'static str;
}

/// Represents discrete Gaussian distribution parameters
#[derive(Debug, Clone, Copy)]
pub struct DiscreteGaussian {
    pub sigma: f64,
    pub bound: f64,
}

impl DistributionParameters for DiscreteGaussian {
    fn type_name(&self) -> &'static str {
        "DiscreteGaussian"
    }
}

/// Represents ternary distribution parameters
#[derive(Debug, Clone, Copy)]
pub struct Ternary {
    pub p: f64,
    pub h: usize,
}

impl DistributionParameters for Ternary {
    fn type_name(&self) -> &'static str {
        "Ternary"
    }
}

/// Represents uniform distribution parameters
#[derive(Debug, Clone, Copy)]
pub struct Uniform;

impl DistributionParameters for Uniform {
    fn type_name(&self) -> &'static str {
        "Uniform"
    }
}

/// Trait for pseudo-random number generators
pub trait PRNG: std::io::Read {
    /// Resets the PRNG to its initial state
    fn reset(&mut self);
}

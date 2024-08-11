use thiserror::Error;

#[derive(Error, Debug)]
pub enum SubRingError {
    #[error("Invalid ring degree: must be a power of 2 greater than {0}")]
    InvalidRingDegree(u64),

    #[error("Invalid modulus for Montgomery reduction: must not be a power of 2 or 0")]
    InvalidModulusForMontgomery,

    #[error("NTT creation failed")]
    NTTCreationFailed,
    #[error("Computation error: {0}")]
    ComputationError(String),
}

#[derive(Error, Debug)]
pub enum RingError {
    #[error("Invalid ring degree {0}: must be a power of 2")]
    InvalidRingDegree(u64),
    #[error("Moduli must be non-empty")]
    EmptyModuli,
    #[error("All moduli must be distinct primes")]
    NonDistinctPrimeModuli,
    #[error("All moduli must be equal to 1 modulo 2N")]
    InvalidModuli,

    // subring errors
    #[error("Subring error: {0}")]
    SubRingError(#[from] SubRingError),
}

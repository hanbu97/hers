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

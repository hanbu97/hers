use thiserror::Error;

#[derive(Error, Debug)]
pub enum SubRingError {
    #[error("Invalid ring degree: must be a power of 2 greater than {0}")]
    InvalidRingDegree(usize),
    #[error("NTT creation failed")]
    NTTCreationFailed,
    #[error("Computation error: {0}")]
    ComputationError(String),
}

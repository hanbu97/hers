use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RingType {
    /// Standard polynomial ring: Z[X]/(X^N + 1)
    /// This is the default ring type, representing polynomials modulo (X^N + 1).
    Standard,

    /// Conjugate-invariant ring: Z[X+X^-1]/(X^2N + 1)
    /// This ring type represents polynomials that are invariant under the conjugation operation.
    ConjugateInvariant,
}

impl ToString for RingType {
    fn to_string(&self) -> String {
        match self {
            RingType::Standard => "Standard".to_string(),
            RingType::ConjugateInvariant => "ConjugateInvariant".to_string(),
        }
    }
}

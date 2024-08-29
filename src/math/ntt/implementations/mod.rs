use crate::math::ntt::{params::NTTTable, traits::NumberTheoreticTransform};

#[derive(Debug, Clone)]
pub struct NTTBase {
    pub degree: usize,       // Degree of the polynomial
    pub modulus: u64,        // Modulus for the ring
    pub m_red_constant: u64, // Montgomery reduction constant
    pub b_red_constant: [u64; 2], // Barrett reduction constants
                             // pub ntt_table: NTTTable,      // Table of NTT constants
}

pub mod standard_ntt;

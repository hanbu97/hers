// Poly struct definition
pub struct Poly {
    pub coeffs: Vec<Vec<u64>>,
}

impl Poly {
    // Constructor for Poly
    pub fn new(levels: usize, n: usize) -> Self {
        Poly {
            coeffs: vec![vec![0; n]; levels],
        }
    }
}

use std::cmp::Ordering;

// Poly struct definition
#[derive(Debug, Clone)]
pub struct Poly {
    pub coeffs: Vec<Vec<u64>>,
}

impl Poly {
    pub fn new(n: usize, level: usize) -> Self {
        let coeffs = vec![vec![0; n]; level + 1];
        Poly { coeffs }
    }

    /// Level returns the current number of moduli minus 1.
    pub fn level(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// Zero sets all coefficients of the target polynomial to 0.
    pub fn zero(&mut self) {
        for coeff_vec in &mut self.coeffs {
            for coeff in coeff_vec.iter_mut() {
                *coeff = 0;
            }
        }
    }

    /// Degree returns the number of coefficients of the polynomial, which equals the degree of the Ring cyclotomic polynomial.
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() {
            0
        } else {
            self.coeffs[0].len()
        }
    }

    // CopyLvl copies the coefficients of 'other' on the target polynomial.
    // This method does nothing if the underlying arrays are the same.
    // Expects the degree of both polynomials to be identical.
    pub fn copy_lvl(&mut self, level: usize, other: &Self) {
        for (self_coeffs, other_coeffs) in self
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter())
            .take(level + 1)
        {
            self_coeffs.clone_from(other_coeffs);
        }
    }

    // Resize resizes the level of the target polynomial to the provided level.
    // If the provided level is larger than the current level, then allocates zero
    // coefficients, otherwise dereferences the coefficients above the provided level.
    pub fn resize(&mut self, level: usize) {
        let degree = self.degree();
        match self.level().cmp(&level) {
            Ordering::Greater => self.coeffs.truncate(level + 1),
            Ordering::Less => {
                let prev_level = self.level();
                self.coeffs
                    .extend((prev_level + 1..=level).map(|_| vec![0; degree]));
            }
            Ordering::Equal => {}
        }
    }

    // Copy copies the coefficients of 'other' on the target polynomial.
    // This method does nothing if the underlying arrays are the same.
    // This method will resize the target polynomial to the level of
    // the input polynomial.
    pub fn copy(&mut self, other: &Self) {
        self.resize(other.level());
        self.copy_lvl(other.level(), other);
    }
}

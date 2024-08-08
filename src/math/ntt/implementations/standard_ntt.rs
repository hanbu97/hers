use crate::math::ntt::params::NTTParams;

use super::*;

pub struct StandardNTT {
    pub base: NTTBase,
}

impl StandardNTT {
    pub fn new(params: NTTParams, ntt_table: NTTTable) -> Self {
        let (n, modulus, _, _, b_red_constant, m_red_constant) = params;
        StandardNTT {
            base: NTTBase {
                n: n as usize,
                modulus,
                m_red_constant,
                b_red_constant,
                ntt_table,
            },
        }
    }
}

impl NumberTheoreticTransform for StandardNTT {
    fn forward(&self, p1: &mut [u64], p2: &[u64]) {
        assert_eq!(p1.len(), self.base.n, "Input length must match NTT size");

        // Placeholder for the actual forward NTT algorithm
        unimplemented!("Forward NTT not yet implemented");
    }

    fn forward_lazy(&self, p1: &mut [u64], p2: &[u64]) {
        unimplemented!("Forward lazy NTT not yet implemented");
    }

    fn backward(&self, p1: &mut [u64], p2: &[u64]) {
        assert_eq!(p1.len(), self.base.n, "Input length must match NTT size");

        // Placeholder for the actual backward NTT algorithm
        unimplemented!("Backward NTT not yet implemented");

        // A real implementation would be similar to forward, but with these changes:
        // 1. Use roots_backward instead of roots_forward
        // 2. Reverse the order of the stages
        // 3. Apply the scaling factor n_inv at the end
    }

    fn backward_lazy(&self, p1: &mut [u64], p2: &[u64]) {
        unimplemented!("Backward lazy NTT not yet implemented");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_roundtrip() {
        // let n = 8;
        // let modulus = 17; // A small prime for testing
        // let ntt = StandardNTT::new(n, modulus);

        // let mut data = vec![1, 2, 3, 4, 5, 6, 7, 8];
        // let original_data = data.clone();

        // // This test will fail until the NTT methods are implemented
        // ntt.forward(&mut data, &[]);
        // ntt.backward(&mut data, &[]);

        // assert_eq!(data, original_data, "NTT roundtrip failed");
    }
}

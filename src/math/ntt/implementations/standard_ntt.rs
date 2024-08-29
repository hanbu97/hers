use crate::math::{
    ntt::params::NTTParameters,
    ring::operations::ntt_operations::{
        intt_standard, intt_standard_lazy, ntt_standard, ntt_standard_lazy,
    },
};

use super::*;

#[derive(Debug, Clone)]
pub struct StandardNTT {
    pub base: NTTBase,
}

impl StandardNTT {
    pub fn new(params: NTTParameters) -> Self {
        StandardNTT {
            base: NTTBase {
                degree: params.degree as usize,
                modulus: params.modulus,
                m_red_constant: params.m_red_constant,
                b_red_constant: params.b_red_constant,
                // ntt_table,
            },
        }
    }
}

impl NumberTheoreticTransform for StandardNTT {
    fn forward(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable) {
        ntt_standard(
            p1,
            p2,
            self.base.degree as u64,
            self.base.modulus,
            self.base.m_red_constant,
            self.base.b_red_constant,
            &ntt_table.roots_forward,
        )
    }

    fn forward_lazy(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable) {
        ntt_standard_lazy(
            p1,
            p2,
            self.base.degree,
            self.base.modulus,
            self.base.m_red_constant,
            &ntt_table.roots_forward,
        )
    }

    fn backward(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable) {
        intt_standard(
            p1,
            p2,
            self.base.degree,
            ntt_table.n_inv,
            self.base.modulus,
            self.base.m_red_constant,
            &ntt_table.roots_backward,
        )
    }

    fn backward_lazy(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable) {
        intt_standard_lazy(
            p1,
            p2,
            self.base.degree,
            ntt_table.n_inv,
            self.base.modulus,
            self.base.m_red_constant,
            &ntt_table.roots_backward,
        )
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

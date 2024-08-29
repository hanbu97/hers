use enum_dispatch::enum_dispatch;

use super::params::NTTTable;

// NTT trait to provide lexibility on what type of NTT is used by the struct Ring.
#[enum_dispatch]
pub trait NumberTheoreticTransform {
    fn forward(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable);
    fn forward_lazy(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable);
    fn backward(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable);
    fn backward_lazy(&self, p1: &[u64], p2: &mut [u64], ntt_table: &NTTTable);
}

use enum_dispatch::enum_dispatch;

use super::params::{NTTParams, NTTTable};

// NTT trait to provide lexibility on what type of NTT is used by the struct Ring.
#[enum_dispatch]
pub trait NumberTheoreticTransform {
    fn forward(&self, p1: &mut [u64], p2: &[u64]);
    fn forward_lazy(&self, p1: &mut [u64], p2: &[u64]);
    fn backward(&self, p1: &mut [u64], p2: &[u64]);
    fn backward_lazy(&self, p1: &mut [u64], p2: &[u64]);
}

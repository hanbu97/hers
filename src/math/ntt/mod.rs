pub mod core;
pub mod implementations;
pub mod params;
pub mod traits;

use enum_dispatch::enum_dispatch;
use traits::NumberTheoreticTransform;

pub use implementations::standard_ntt::StandardNTT;

use self::params::{NTTParams, NTTTable};

// In ntt/implementations/mod.rs
#[enum_dispatch(NumberTheoreticTransform)]
pub enum NTTImplementations {
    StandardNTT,
    // OptimizedNTT,
    // Add other NTT implementations here
}

impl NTTImplementations {
    pub fn new_standard(params: NTTParams, ntt_table: NTTTable) -> Self {
        NTTImplementations::StandardNTT(StandardNTT::new(params, ntt_table))
    }
}

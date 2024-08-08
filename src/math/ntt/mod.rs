pub mod implementations;
pub mod params;
pub mod traits;

use enum_dispatch::enum_dispatch;
use traits::NumberTheoreticTransform;

pub use implementations::standard_ntt::StandardNTT;

// In ntt/implementations/mod.rs
#[enum_dispatch(NumberTheoreticTransform)]
pub enum NTTImplementations {
    StandardNTT,
    // OptimizedNTT,
    // Add other NTT implementations here
}

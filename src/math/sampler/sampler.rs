use enum_dispatch::enum_dispatch;
use rand_core::RngCore;

use crate::math::ring::{polynomial::Poly, Ring};

use super::{
    gaussian_sampler::{DiscreteGaussian, GaussianSampler},
    ternary_sampler::{Ternary, TernarySampler},
    uniform_sampler::UniformSampler,
};

#[enum_dispatch(SamplerEnum)]
pub trait Sampler {
    type Rng: RngCore + Clone;

    fn read(&mut self, pol: &mut Poly);
    fn read_new(&mut self) -> Poly;
    fn read_and_add(&mut self, pol: &mut Poly);
    fn at_level(&self, level: usize) -> SamplerEnum<Self::Rng>;
}

pub enum SamplerEnum<R: RngCore + Clone> {
    Gaussian(GaussianSampler<R>),
    Ternary(TernarySampler<R>),
    Uniform(UniformSampler<R>),
}

pub enum DistributionParameters {
    DiscreteGaussian { sigma: f64, bound: f64 },
    Ternary { p: f64, h: usize },
    Uniform,
}

pub fn new_sampler<R: RngCore + Clone + 'static>(
    prng: R,
    base_ring: Ring,
    params: DistributionParameters,
    montgomery: bool,
) -> Result<SamplerEnum<R>, String> {
    match params {
        DistributionParameters::DiscreteGaussian { sigma, bound } => {
            Ok(SamplerEnum::Gaussian(GaussianSampler::new(
                prng,
                base_ring,
                DiscreteGaussian { sigma, bound },
                montgomery,
            )))
        }
        DistributionParameters::Ternary { p, h } => {
            let ternary = if p != 0.0 {
                Ternary::Probability(p)
            } else {
                Ternary::HammingWeight(h)
            };
            Ok(SamplerEnum::Ternary(TernarySampler::new(
                prng, base_ring, ternary, montgomery,
            )?))
        }
        DistributionParameters::Uniform => {
            Ok(SamplerEnum::Uniform(UniformSampler::new(prng, base_ring)))
        }
    }
}

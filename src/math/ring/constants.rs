/// Constants related to Galois theory and ring operations in homomorphic encryption.

/// The generator for Galois automorphisms.
///
/// This constant is used in the context of ring automorphisms in homomorphic encryption schemes.
///
/// Properties:
/// - It is an integer of order N/2 modulo M, where N is the ring degree and M is the coefficient modulus.
/// - It spans Z_M (the integers modulo M) together with -1.
/// - The j-th ring automorphism takes the root zeta to zeta^(5^j).
///
/// Choice of 5:
/// - 5 is the smallest odd prime that satisfies the required properties for most practical ring degrees.
/// - It simplifies computations and is widely used in many homomorphic encryption implementations.
/// - For a ring degree N, 5 generates the multiplicative subgroup of order N/2 modulo 2N.
///
/// Note: While 5 is a common and efficient choice, other values could theoretically be used
/// if they satisfy the necessary mathematical properties.
pub const GALOIS_GEN: u64 = 5;

/// The minimum ring degree required for safely performing loop-unrolled operations.
///
/// This constant defines the smallest ring size that allows for efficient loop unrolling
/// and SIMD (Single Instruction, Multiple Data) optimizations.
///
/// Reasons for choosing 8:
/// 1. It's the first power of 2 large enough to effectively utilize modern CPU's SIMD instructions.
/// 2. For degrees smaller than 8, loop unrolling might not provide significant performance benefits.
/// 3. 8 is generally friendly for memory alignment and cache line sizes on most architectures.
///
/// Usage:
/// - Serves as a threshold to determine when to apply certain optimization techniques.
/// - Ensures that there's enough data to fully utilize parallel processing capabilities.
/// - Prevents over-optimization for very small data sets, which could lead to performance degradation.
///
/// Note: This value might need adjustment based on specific hardware architectures or if
/// future CPUs have different optimal sizes for SIMD operations.
pub const MINIMUM_RING_DEGREE_FOR_LOOP_UNROLLED_OPS: usize = 8;

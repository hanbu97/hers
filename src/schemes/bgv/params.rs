// use std::sync::Arc;

// /// Represents the parameters for the BGV (Brakerski-Gentry-Vaikuntanathan) homomorphic encryption scheme.
// /// This structure encapsulates all necessary parameters for key generation, encryption, decryption,
// /// and homomorphic operations in the BGV scheme.
// pub struct BGVParameters {
//     /// Core parameters for the Ring-Learning With Errors (RLWE) problem.
//     /// These form the cryptographic basis for the BGV scheme's security.
//     rlwe_params: RLWEParameters,

//     /// Polynomial ring used for multiplication operations in the ciphertext space.
//     /// This ring is defined over Z_q[X]/(X^N + 1), where q is the ciphertext modulus
//     /// and N is the ring degree.
//     /// It's used to optimize homomorphic multiplication operations.
//     ring_q_mul: Arc<Ring>,

//     /// Polynomial ring representing the plaintext space.
//     /// This ring is defined over Z_t[X]/(X^N + 1), where t is the plaintext modulus
//     /// and N is the ring degree.
//     /// It's used for encoding and decoding messages, and for plaintext operations.
//     ring_t: Arc<Ring>,
// }

// impl BGVParameters {
//     /// Creates a new instance of BGVParameters.
//     ///
//     /// # Arguments
//     ///
//     /// * `rlwe_params` - The RLWE parameters.
//     /// * `ring_q_mul` - The ring for ciphertext multiplication.
//     /// * `ring_t` - The ring for plaintext operations.
//     ///
//     /// # Returns
//     ///
//     /// A new `BGVParameters` instance.
//     pub fn new(rlwe_params: RLWEParameters, ring_q_mul: Ring, ring_t: Ring) -> Self {
//         Self {
//             rlwe_params,
//             ring_q_mul: Arc::new(ring_q_mul),
//             ring_t: Arc::new(ring_t),
//         }
//     }

//     /// Returns a reference to the RLWE parameters.
//     pub fn rlwe_params(&self) -> &RLWEParameters {
//         &self.rlwe_params
//     }

//     /// Returns a reference to the ciphertext multiplication ring.
//     pub fn ring_q_mul(&self) -> &Arc<Ring> {
//         &self.ring_q_mul
//     }

//     /// Returns a reference to the plaintext ring.
//     pub fn ring_t(&self) -> &Arc<Ring> {
//         &self.ring_t
//     }
// }

use num_bigint::BigUint;
use num_traits::ToPrimitive;

/// Computes the constant for the BRed algorithm.
/// Returns [((2^128)/q)/(2^64), (2^128)/q mod 2^64].
pub fn compute_barrett_constants(p: u64) -> [u64; 2] {
    let barrett = ((BigUint::from(1u64) << 128usize) / BigUint::from(p))
        .to_u128()
        .unwrap();
    [(barrett >> 64) as u64, barrett as u64]
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_barrett_constants() {
        let test_cases = [
            (3, 6148914691236517205, 6148914691236517205),
            (17, 1085102592571150095, 1085102592571150095),
            (257, 71777214294589695, 71777214294589695),
            (65537, 281470681808895, 281470681808895),
            (4294967295, 4294967297, 4294967297),
            (4294967296, 4294967296, 0),
            (9223372036854775807, 2, 4),
            (18446744073709551615, 1, 1),
        ];

        for &(q, expected_hi, expected_lo) in &test_cases {
            let constants = compute_barrett_constants(q);
            println!("q: {}", q);
            println!("((2^128)/q)/(2^64): {}", constants[0]);
            println!("(2^128)/q mod 2^64: {}\n", constants[1]);

            assert_eq!(
                constants[0], expected_hi,
                "Mismatch in high 64 bits for q = {}",
                q
            );
            assert_eq!(
                constants[1], expected_lo,
                "Mismatch in low 64 bits for q = {}",
                q
            );
        }
    }
}

/// Computes the constant m_red_constant = (q^-1) mod 2^64 required for MRed.
pub fn compute_montgomery_constant(q: u64) -> u64 {
    let mut mredconstant = 1u64;
    let mut q_temp = q;
    for _ in 0..63 {
        mredconstant = mredconstant.wrapping_mul(q_temp);
        q_temp = q_temp.wrapping_mul(q_temp);
    }
    mredconstant
}

#[test]
fn test_compute_montgomery_constant() {
    let test_cases = [
        (3, 12297829382473034411),
        (17, 17361641481138401521),
        (257, 18374966859414961921),
        (65537, 18446462603027742721),
        (4294967295, 18446744069414584319),           // 2^32 - 1
        (4294967296, 0),                              // 2^32
        (9223372036854775807, 9223372036854775807),   // 2^63 - 1
        (18446744073709551615, 18446744073709551615), // 2^64 - 1
    ];

    for &(q, expected_m_red_constant) in &test_cases {
        let m_red_constant = compute_montgomery_constant(q);

        // println!("q: {}", q);
        // println!("m_red_constant: {}", m_red_constant);
        // println!("expected: {}\n", expected_m_red_constant);

        assert_eq!(
            m_red_constant, expected_m_red_constant,
            "Mismatch for q = {}: expected {}, but got {}",
            q, expected_m_red_constant, m_red_constant
        );
    }
}

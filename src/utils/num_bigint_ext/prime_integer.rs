use num_traits::FromPrimitive;
use num_traits::One;
use rug::integer::IntegerExt64;
use rug::rand::RandState;

use rug::Integer;

use super::constants_integer::{BIG_1, BIG_2, BIG_3, BIG_64, PRIMES_A, PRIMES_B, PRIME_BIT_MASK};

pub fn probably_prime_lucas(n: &Integer) -> bool {
    if n.is_zero() || n.is_one() {
        return false;
    }

    // Two is the only even prime.
    if n.to_u64() == Some(2) {
        return false;
    }

    let mut p = 3u64;
    let n_int = n.clone();

    loop {
        if p > 10000 {
            // This is widely believed to be impossible.
            // If we get a report, we'll want the exact number n.
            panic!("internal error: cannot find (D/n) = -1 for {:?}", n)
        }

        let d_int = Integer::from_u64(p * p - 4).unwrap();
        let j = Integer::jacobi(&d_int, &n_int);

        if j == -1 {
            break;
        }
        if j == 0 {
            // d = p²-4 = (p-2)(p+2).
            // If (d/n) == 0 then d shares a prime factor with n.
            // Since the loop proceeds in increasing p and starts with p-2==1,
            // the shared prime factor must be p+2.
            // If p+2 == n, then n is prime; otherwise p+2 is a proper factor of n.
            return n_int.to_i64() == Some(p as i64 + 2);
        }
        if p == 40 {
            // We'll never find (d/n) = -1 if n is a square.
            // If n is a non-square we expect to find a d in just a few attempts on average.
            // After 40 attempts, take a moment to check if n is indeed a square.
            let t1 = n.clone().sqrt();
            if Integer::from(&t1 * &t1) == *n {
                return false;
            }
        }

        p += 1;
    }

    let mut s: Integer = Integer::from(n + 1);
    let r = s.find_one(0).unwrap();
    s >>= r;

    let nm2 = Integer::from(n - &*BIG_2); // n - 2

    let mut vk = BIG_2.clone();
    let mut vk1 = Integer::from_u64(p).unwrap();

    for i in (0..s.significant_bits_64()).rev() {
        if s.get_bit_64(i) {
            // k' = 2k+1
            // V(k') = V(2k+1) = V(k) V(k+1) - P
            // let t1 = Integer::from((&vk * &vk1) + n) - p.into();
            let t1: Integer = (&vk * &vk1 - &p.into()).into();
            vk = t1 % n;
            // V(k'+1) = V(2k+2) = V(k+1)² - 2
            let t1 = Integer::from(&vk1 * &vk1) + &nm2;
            vk1 = t1 % n;
        } else {
            // k' = 2k
            // V(k'+1) = V(2k+1) = V(k) V(k+1) - P
            let t1 = Integer::from((&vk * &vk1) + n) - p;
            vk1 = t1 % n;
            // V(k') = V(2k) = V(k)² - 2
            let t1 = Integer::from((&vk * &vk) + &nm2);
            vk = t1 % n;
        }
    }

    // Now k=s, so vk = V(s). Check V(s) ≡ ±2 (mod n).
    if vk.to_u64() == Some(2) || vk == nm2 {
        // Check U(s) ≡ 0.
        // As suggested by Jacobsen, apply Crandall and Pomerance equation 3.13:
        //
        //	U(k) = D⁻¹ (2 V(k+1) - P V(k))
        //
        // Since we are checking for U(k) == 0 it suffices to check 2 V(k+1) == P V(k) mod n,
        // or P V(k) - 2 V(k+1) == 0 mod n.
        let mut t1 = Integer::from(&vk * p);
        let mut t2 = Integer::from(&vk1 << 1);

        if t1 < t2 {
            core::mem::swap(&mut t1, &mut t2);
        }

        t1 -= t2;

        if (t1 % n).is_zero() {
            return true;
        }
    }

    // Check V(2^t s) ≡ 0 mod n for some 0 ≤ t < r-1.
    for _ in 0..r - 1 {
        if vk.is_zero() {
            return true;
        }

        // Optimization: V(k) = 2 is a fixed point for V(k') = V(k)² - 2,
        // so if V(k) = 2, we can stop: we will never find a future V(k) == 0.
        if vk.to_u64() == Some(2) {
            return false;
        }

        // k' = 2k
        // V(k') = V(2k) = V(k)² - 2
        let t1 = Integer::from((&vk * &vk) - &*BIG_2);
        vk = t1 % n;
    }

    false
}

pub fn probably_prime_miller_rabin(n: &Integer, reps: usize, force2: bool) -> bool {
    let nm1: Integer = (n - &*BIG_1).into();
    let k = nm1.find_one(0).unwrap();
    let q: Integer = (&nm1 >> k).into();

    let nm3: Integer = (n - &*BIG_3).into();

    let mut rng = RandState::new();

    'nextrandom: for i in 0..reps {
        let x = if i == reps - 1 && force2 {
            Integer::from(2)
        } else {
            let rand_bits = nm3.significant_bits();
            Integer::from(rng.bits(rand_bits as u32)) % &nm3 + 2
        };

        let mut y = x.pow_mod(&q, n).unwrap();
        if y == 1 || y == nm1 {
            continue;
        }

        for _ in 1..k {
            y = y.pow_mod(&Integer::from(2), n).unwrap();
            if y == nm1 {
                continue 'nextrandom;
            }
            if y == 1 {
                return false;
            }
        }
        return false;
    }

    true
}

pub fn probably_prime(x: &Integer, n: usize) -> bool {
    if x.is_zero() {
        return false;
    }

    if x < &*BIG_64 {
        return (PRIME_BIT_MASK & (1u64 << x.to_u32().unwrap())) != 0;
    }

    if x.is_even() {
        return false;
    }

    let r_a = Integer::from(x % &PRIMES_A);
    let r_b = Integer::from(x % &PRIMES_B);

    if r_a.is_divisible(&Integer::from(3u32))
        || r_a.is_divisible(&Integer::from(5u32))
        || r_a.is_divisible(&Integer::from(7u32))
        || r_a.is_divisible(&Integer::from(11u32))
        || r_a.is_divisible(&Integer::from(13u32))
        || r_a.is_divisible(&Integer::from(17u32))
        || r_a.is_divisible(&Integer::from(19u32))
        || r_a.is_divisible(&Integer::from(23u32))
        || r_a.is_divisible(&Integer::from(37u32))
        || r_b.is_divisible(&Integer::from(29u32))
        || r_b.is_divisible(&Integer::from(31u32))
        || r_b.is_divisible(&Integer::from(41u32))
        || r_b.is_divisible(&Integer::from(43u32))
        || r_b.is_divisible(&Integer::from(47u32))
        || r_b.is_divisible(&Integer::from(53u32))
    {
        return false;
    }

    probably_prime_miller_rabin(x, n + 1, true) && probably_prime_lucas(x)
}

#[cfg(test)]
mod tests {
    use crate::utils::prime::PrimeChecking;

    use super::*;

    lazy_static::lazy_static! {
        pub static ref PRIMES: Vec<&'static str> = vec![
        "2",
        "3",
        "5",
        "7",
        "11",

        "13756265695458089029",
        "13496181268022124907",
        "10953742525620032441",
        "17908251027575790097",

        // https://golang.org/issue/638
        "18699199384836356663",

        "98920366548084643601728869055592650835572950932266967461790948584315647051443",
        "94560208308847015747498523884063394671606671904944666360068158221458669711639",

        // http://primes.utm.edu/lists/small/small3.html
        "449417999055441493994709297093108513015373787049558499205492347871729927573118262811508386655998299074566974373711472560655026288668094291699357843464363003144674940345912431129144354948751003607115263071543163",
        "230975859993204150666423538988557839555560243929065415434980904258310530753006723857139742334640122533598517597674807096648905501653461687601339782814316124971547968912893214002992086353183070342498989426570593",
        "5521712099665906221540423207019333379125265462121169655563495403888449493493629943498064604536961775110765377745550377067893607246020694972959780839151452457728855382113555867743022746090187341871655890805971735385789993",
        "203956878356401977405765866929034577280193993314348263094772646453283062722701277632936616063144088173312372882677123879538709400158306567338328279154499698366071906766440037074217117805690872792848149112022286332144876183376326512083574821647933992961249917319836219304274280243803104015000563790123",
        // ECC primes: http://tools.ietf.org/html/draft-ladd-safecurves-02
        "3618502788666131106986593281521497120414687020801267626233049500247285301239",                                                                                  // Curve1174: 2^251-9
        "57896044618658097711785492504343953926634992332820282019728792003956564819949",                                                                                 // Curve25519: 2^255-19
        "9850501549098619803069760025035903451269934817616361666987073351061430442874302652853566563721228910201656997576599",                                           // E-382: 2^382-105
        "42307582002575910332922579714097346549017899709713998034217522897561970639123926132812109468141778230245837569601494931472367",                                 // Curve41417: 2^414-17
        "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151", // E-521: 2^521-1
        ];

        static ref COMPOSITES: Vec<&'static str> = vec![
            "0",
            "1",

            "21284175091214687912771199898307297748211672914763848041968395774954376176754",
            "6084766654921918907427900243509372380954290099172559290432744450051395395951",
            "84594350493221918389213352992032324280367711247940675652888030554255915464401",
            "82793403787388584738507275144194252681",

            // Arnault, "Rabin-Miller Primality Test: Composite Numbers Which Pass It",
            // Mathematics of Computation, 64(209) (January 1995), pp. 335-361.
            "1195068768795265792518361315725116351898245581", // strong pseudoprime to prime bases 2 through 29
            // strong pseudoprime to all prime bases up to 200
            "8038374574536394912570796143419421081388376882875581458374889175222974273765333652186502336163960045457915042023603208766569966760987284043965408232928738791850869166857328267761771029389697739470167082304286871099974399765441448453411558724506334092790222752962294149842306881685404326457534018329786111298960644845216191652872597534901",

            // Extra-strong Lucas pseudoprimes. https://oeis.org/A217719
            "989",
            "3239",
            "5777",
            "10877",
            "27971",
            "29681",
            "30739",
            "31631",
            "39059",
            "72389",
            "73919",
            "75077",
            "100127",
            "113573",
            "125249",
            "137549",
            "137801",
            "153931",
            "155819",
            "161027",
            "162133",
            "189419",
            "218321",
            "231703",
            "249331",
            "370229",
            "429479",
            "430127",
            "459191",
            "473891",
            "480689",
            "600059",
            "621781",
            "632249",
            "635627",

            "3673744903",
            "3281593591",
            "2385076987",
            "2738053141",
            "2009621503",
            "1502682721",
            "255866131",
            "117987841",
            "587861",

            "6368689",
            "8725753",
            "80579735209",
            "105919633",
        ];

        // Test Cases from #51
        static ref ISSUE_51: Vec<&'static str> = vec![
            "1579751",
            "1884791",
            "3818929",
            "4080359",
            "4145951",
        ];
    }

    #[test]
    fn test_bit_set() {
        let v = vec![0b10101001u32];
        let num = Integer::from_digits(&v, rug::integer::Order::Lsf);
        assert!(num.get_bit(0));
        assert!(!num.get_bit(1));
        assert!(!num.get_bit(2));
        assert!(num.get_bit(3));
        assert!(!num.get_bit(4));
        assert!(num.get_bit(5));
        assert!(!num.get_bit(6));
        assert!(num.get_bit(7));
    }

    #[test]
    fn test_next_prime_basics() {
        let primes1 = (0..2048u32)
            .map(|i| Integer::from(i).next_prime())
            .collect::<Vec<_>>();
        let primes2 = (0..2048u32)
            .map(|i| {
                let i = Integer::from(i);
                let p = i.clone().next_prime();
                assert!(p > i);
                p
            })
            .collect::<Vec<_>>();

        for (p1, p2) in primes1.iter().zip(&primes2) {
            assert_eq!(p1, p2);
            assert!(p1.is_prime());
        }
    }

    #[test]
    fn test_next_prime_bug_44() {
        let i = Integer::from(1032989);
        let next = i.next_prime();
        assert_eq!(Integer::from(1033001), next);
    }

    macro_rules! test_pseudo_primes {
        ($name:ident, $cond:expr, $want:expr) => {
            #[test]
            fn $name() {
                let mut i = 3;
                let mut want = $want;
                while i < 100000 {
                    let n = Integer::from(i);
                    let pseudo = $cond(&n);
                    if pseudo && (want.is_empty() || i != want[0]) {
                        panic!("cond({}) = true, want false", i);
                    } else if !pseudo && !want.is_empty() && i == want[0] {
                        panic!("cond({}) = false, want true", i);
                    }
                    if !want.is_empty() && i == want[0] {
                        want = want[1..].to_vec();
                    }
                    i += 2;
                }

                if !want.is_empty() {
                    panic!("forgot to test: {:?}", want);
                }
            }
        };
    }

    test_pseudo_primes!(
        test_probably_prime_miller_rabin,
        |n| probably_prime_miller_rabin(n, 1, true) && !probably_prime_lucas(n),
        vec![
            2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665, 80581,
            85489, 88357, 90751,
        ]
    );

    test_pseudo_primes!(
        test_probably_prime_lucas,
        |n| probably_prime_lucas(n) && !probably_prime_miller_rabin(n, 1, true),
        vec![989, 3239, 5777, 10877, 27971, 29681, 30739, 31631, 39059, 72389, 73919, 75077,]
    );
}

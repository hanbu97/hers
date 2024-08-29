pub mod core;
pub mod implementations;
pub mod params;
pub mod traits;

use enum_dispatch::enum_dispatch;
use traits::NumberTheoreticTransform;

pub use implementations::standard_ntt::StandardNTT;

use self::params::{NTTParameters, NTTTable};

// In ntt/implementations/mod.rs
#[enum_dispatch(NumberTheoreticTransform)]
#[derive(Debug, Clone)]
pub enum NTTImplementations {
    StandardNTT,
    // OptimizedNTT,
    // Add other NTT implementations here
}

impl NTTImplementations {
    pub fn new_standard(params: NTTParameters, ntt_table: NTTTable) -> Self {
        NTTImplementations::StandardNTT(StandardNTT::new(params, ntt_table))
    }
}

#[cfg(test)]
mod test {
    use crate::math::ring::{polynomial::Poly, Ring};

    // use super::*;

    struct TestCase {
        pub degree: u64,
        pub qis: Vec<u64>,
        pub poly: Poly,
        pub poly_ntt: Poly,
    }

    lazy_static::lazy_static! {
        static ref TEST_VECTOR: Vec<TestCase> = vec![
            TestCase {
                degree: 16,
                qis: vec![576460752303439873, 576460752303702017],
                poly: Poly {
                    coeffs: vec![
                        vec![29335002291498019, 74733314878908829, 345757914625392883, 424592696763883150, 305098757618029540, 315880659253740539, 566291353020324899, 381879490285643315, 34642655966258078, 436368737741273744, 422320479487058982, 251503834452711492, 379754966293786644, 266993967580766257, 265441209649369663, 479048496297441983],
                        vec![229005636957624603, 39991394218169426, 168047666046761487, 148360907414915405, 73259769245767872, 16981974422312794, 496977853225992141, 166066041724987771, 264052080009592093, 298274702686123828, 35777507392976624, 357559017452722394, 314515717429384298, 162821044855043426, 109977030677147798, 81303063671114932],
                    ]
                },
                poly_ntt: Poly {
                    coeffs: vec![
                        vec![478709994917861263, 384523361984839039, 85280178929118517, 97236771105538581, 405398446277957930, 212032954159995430, 422470404160315474, 554803939008707088, 548834797847219388, 77555291080479046, 395019082584063204, 199181437220481637, 117237287301343342, 288680759037675256, 399758453229973389, 414322896245918704],
                        vec![48052203194603178, 560437377430510021, 51924270083317129, 254030332439706305, 520426933791709415, 443676955646482348, 405741025864202685, 70579349438930370, 187051495725458514, 84142641467084820, 194371127241444851, 191269223870154261, 109044160236534164, 304031719544775780, 243823945337031160, 571948182313750664],
                    ]
                },
            },
        ];
    }

    #[test]
    fn test_ntt() {
        for test_case in TEST_VECTOR.iter() {
            let ring =
                Ring::new(test_case.degree, test_case.qis.clone()).expect("Failed to create ring");

            println!(
                "degree: {}, limbs: {}",
                ring.degree(),
                ring.moduli_chain_length()
            );

            let mut x = ring.new_poly();
            let mut y = ring.new_poly();
            let mut z = ring.new_poly();

            x.copy(&test_case.poly);
            y.copy(&test_case.poly_ntt);

            // ring.ntt(&x, &mut z);

            println!("subrings len: {}", ring.sub_rings.len());
            println!(
                "ring ntt roots: {:?}",
                ring.sub_rings[0].ntt_table.roots_forward
            );

            println!(
                "ring ntt roots: {:?}",
                ring.sub_rings[1].ntt_table.roots_forward
            );

            // println!("x: {:?}", x);
            // println!("y: {:?}", y);
            // println!("z: {:?}", z);
        }
    }
}

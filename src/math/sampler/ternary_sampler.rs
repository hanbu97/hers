use itertools::Itertools;
use rand_core::RngCore;
use std::f64::consts::E;

use crate::math::ring::{polynomial::Poly, Ring};

const TERNARY_SAMPLER_PRECISION: usize = 56;

pub struct TernarySampler<R: RngCore + Clone> {
    base_ring: Ring,
    prng: R,
    matrix_proba: [[u8; TERNARY_SAMPLER_PRECISION - 1]; 2],
    matrix_values: Vec<[u64; 3]>,
    inv_density: f64,
    hw: usize,
    montgomery: bool,
}

pub enum Ternary {
    Probability(f64),
    HammingWeight(usize),
}

impl<R: RngCore + Clone> TernarySampler<R> {
    pub fn new(prng: R, base_ring: Ring, x: Ternary, montgomery: bool) -> Result<Self, String> {
        println!("Creating new TernarySampler");
        println!("montgomery: {}", montgomery);
        let mut ts = TernarySampler {
            base_ring,
            prng,
            matrix_proba: [[0; TERNARY_SAMPLER_PRECISION - 1]; 2],
            matrix_values: Vec::new(),
            inv_density: 0.0,
            hw: 0,
            montgomery,
        };

        ts.initialize_matrix();

        match x {
            Ternary::Probability(p) if p != 0.0 => {
                ts.inv_density = 1.0 - p;

                println!("Probability: {}, inv_density: {}", p, ts.inv_density);

                if ts.inv_density != 0.5 {
                    ts.compute_matrix_ternary(ts.inv_density);
                }
            }
            Ternary::HammingWeight(h) if h != 0 => {
                ts.hw = h;
            }
            _ => {
                return Err(
                    "Invalid TernarySampler: exactly one of (H, P) should be > 0".to_string(),
                )
            }
        }

        Ok(ts)
    }

    pub fn at_level(&self, level: usize) -> Self {
        TernarySampler {
            base_ring: self.base_ring.at_level(level),
            prng: self.prng.clone(),
            matrix_proba: self.matrix_proba,
            matrix_values: self.matrix_values.clone(),
            inv_density: self.inv_density,
            hw: self.hw,
            montgomery: self.montgomery,
        }
    }

    pub fn read(&mut self, pol: &mut Poly) {
        if self.hw != 0 {
            self.sample_sparse(pol, |_, b, _| b);
        } else {
            self.sample_proba(pol, |_, b, _| b);
        }
    }

    pub fn read_new(&mut self) -> Poly {
        let mut pol = self.base_ring.new_poly();
        self.read(&mut pol);
        pol
    }

    pub fn read_and_add(&mut self, pol: &mut Poly) {
        if self.hw != 0 {
            self.sample_sparse(pol, |a, b, c| (a + b) % c);
        } else {
            self.sample_proba(pol, |a, b, c| (a + b) % c);
        }
    }

    fn initialize_matrix(&mut self) {
        self.matrix_values = self
            .base_ring
            .moduli_chain()
            .iter()
            .enumerate()
            .map(|(i, &modulus)| {
                let mut values = [0u64; 3];
                values[0] = 0;
                if self.montgomery {
                    let mut one_poly = self.base_ring.new_poly();
                    let mut minus_one_poly = self.base_ring.new_poly();
                    one_poly.coeffs[i][0] = 1;
                    minus_one_poly.coeffs[i][0] = modulus - 1;

                    let mut one_m_form = self.base_ring.new_poly();
                    let mut minus_one_m_form = self.base_ring.new_poly();

                    self.base_ring.m_form(&one_poly, &mut one_m_form);
                    self.base_ring
                        .m_form(&minus_one_poly, &mut minus_one_m_form);

                    values[1] = one_m_form.coeffs[i][0];
                    values[2] = minus_one_m_form.coeffs[i][0];
                } else {
                    values[1] = 1;
                    values[2] = modulus - 1;
                }
                values
            })
            .collect();
    }

    fn compute_matrix_ternary(&mut self, p: f64) {
        let g1 = p * E.powi(TERNARY_SAMPLER_PRECISION as i32);
        let x1 = g1 as u64;
        for j in 0..TERNARY_SAMPLER_PRECISION - 1 {
            self.matrix_proba[0][j] = ((x1 >> (TERNARY_SAMPLER_PRECISION - j - 1)) & 1) as u8;
        }

        let g2 = (1.0 - p) * E.powi(TERNARY_SAMPLER_PRECISION as i32);
        let x2 = g2 as u64;
        for j in 0..TERNARY_SAMPLER_PRECISION - 1 {
            self.matrix_proba[1][j] = ((x2 >> (TERNARY_SAMPLER_PRECISION - j - 1)) & 1) as u8;
        }
    }

    fn sample_proba<F>(&mut self, pol: &mut Poly, f: F)
    where
        F: Fn(u64, u64, u64) -> u64,
    {
        println!("Entering sample_proba");
        println!("inv_density: {}", self.inv_density);

        let n = self.base_ring.degree() as usize;
        let moduli = self.base_ring.moduli_chain();

        if self.inv_density == 0.5 {
            let mut random_bytes_coeffs = vec![0u8; n >> 3];
            let mut random_bytes_sign = vec![0u8; n >> 3];
            self.prng.fill_bytes(&mut random_bytes_coeffs);
            self.prng.fill_bytes(&mut random_bytes_sign);

            for i in 0..n {
                let coeff = (random_bytes_coeffs[i >> 3] >> (i & 7)) & 1;
                let sign = (random_bytes_sign[i >> 3] >> (i & 7)) & 1;
                let index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1);

                if i % 100 == 0 {
                    println!(
                        "i: {}, coeff: {}, sign: {}, index: {}",
                        i, coeff, sign, index
                    );
                }

                for (j, &qi) in moduli.iter().enumerate() {
                    let old_value = pol.coeffs[j][i];
                    let lut_value = self.matrix_values[j][index as usize];
                    let new_value = f(old_value, lut_value, qi);
                    pol.coeffs[j][i] = new_value;
                }
            }
        } else {
            let mut random_bytes = vec![0u8; n];
            self.prng.fill_bytes(&mut random_bytes);

            println!(
                "random_bytes: {:?}",
                &random_bytes[random_bytes.len() - 20..]
            );

            let mut pointer = 0u8;
            let mut byte_pointer = 0;

            for i in 0..n {
                let (coeff, sign, new_bytes, new_pointer, new_byte_pointer) =
                    self.kysampling(&mut random_bytes, pointer, byte_pointer, n);
                random_bytes = new_bytes;
                pointer = new_pointer;
                byte_pointer = new_byte_pointer;

                let index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1);

                if i % 100 == 0 {
                    println!(
                        "i: {}, coeff: {}, sign: {}, index: {}",
                        i, coeff, sign, index
                    );
                }

                for (j, &qi) in moduli.iter().enumerate() {
                    pol.coeffs[j][i] =
                        f(pol.coeffs[j][i], self.matrix_values[j][index as usize], qi);
                }
            }
        }
    }

    fn sample_sparse<F>(&mut self, pol: &mut Poly, f: F)
    where
        F: Fn(u64, u64, u64) -> u64,
    {
        let n = self.base_ring.degree();
        let hw = std::cmp::min(self.hw, n as usize);
        let moduli = self.base_ring.moduli_chain();

        let mut index: Vec<usize> = (0..n as usize).collect();
        let mut random_bytes = vec![0u8; (hw as f64 / 8.0).ceil() as usize];
        self.prng.fill_bytes(&mut random_bytes);
        let mut pointer = 0u8;

        for i in 0..hw {
            let mask = (1u64 << (64 - (n as usize - i).leading_zeros())) - 1;
            let mut j = self.rand_int32(mask) as usize;
            while j >= n as usize - i {
                j = self.rand_int32(mask) as usize;
            }

            let coeff = (random_bytes[0] >> (i & 7)) & 1;
            let idx_j = index[j];

            for (k, &qi) in moduli.iter().enumerate() {
                pol.coeffs[k][idx_j] = f(
                    pol.coeffs[k][idx_j],
                    self.matrix_values[k][(coeff + 1) as usize],
                    qi,
                );
            }

            index.swap_remove(j);

            pointer += 1;
            if pointer == 8 {
                random_bytes.remove(0);
                pointer = 0;
            }
        }

        for &i in &index {
            for coeffs in pol.coeffs.iter_mut() {
                coeffs[i] = 0;
            }
        }
    }

    fn kysampling(
        &mut self,
        random_bytes: &mut Vec<u8>,
        pointer: u8,
        byte_pointer: usize,
        byte_length: usize,
    ) -> (u64, u64, Vec<u8>, u8, usize) {
        let mut d = 0;
        let mut col = 0;
        let col_len = self.matrix_proba.len();
        let mut current_pointer = pointer;
        let mut current_byte_pointer = byte_pointer;

        loop {
            for i in current_pointer..8 {
                d = (d << 1) + 1 - ((random_bytes[current_byte_pointer] >> i) & 1) as i32;

                if d > col_len as i32 - 1 {
                    self.prng.fill_bytes(random_bytes);
                    return self.kysampling(random_bytes, 0, 0, byte_length);
                }

                for row in (0..col_len).rev() {
                    d -= self.matrix_proba[row][col] as i32;

                    if d == -1 {
                        let sign = if i == 7 {
                            current_pointer = 0;
                            current_byte_pointer += 1;
                            if current_byte_pointer >= byte_length {
                                current_byte_pointer = 0;
                                self.prng.fill_bytes(random_bytes);
                            }
                            random_bytes[current_byte_pointer] & 1
                        } else {
                            current_pointer = i;
                            (random_bytes[current_byte_pointer] >> (i + 1)) & 1
                        };

                        return (
                            row as u64,
                            sign as u64,
                            random_bytes.to_vec(),
                            current_pointer + 1,
                            current_byte_pointer,
                        );
                    }
                }

                col += 1;
            }

            current_pointer = 0;
            current_byte_pointer += 1;

            if current_byte_pointer >= byte_length {
                current_byte_pointer = 0;
                self.prng.fill_bytes(random_bytes);
            }
        }
    }

    fn rand_int32(&mut self, mask: u64) -> u64 {
        let mut random_bytes = [0u8; 4];
        self.prng.fill_bytes(&mut random_bytes);
        u32::from_be_bytes(random_bytes) as u64 & mask
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::rand::fixed_rand::FixedRandom;

    use super::*;
    // use rand::SeedableRng;
    // use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_ternary_sampler_consistency() {
        let random_numbers_prob: Vec<u64> = vec![14979877486494237012];

        let log_n = 10;
        let qi: Vec<u64> = vec![
            2305843009137934337,
            // 2305843009132953601,
            // 2305843009131642881,
            // 2305843009127448577,
            // 2305843009117224961,
            // 2305843009114341377,
            // 2305843009110409217,
            // 2305843009097564161,
            // 2305843009090748417,
            // 2305843009090486273,
            // 2305843009081311233,
            // 2305843009071087617,
            // 2305843009069514753,
            // 2305843009063223297,
        ];

        let ring = Ring::new(1 << log_n, qi).unwrap();

        let rands: Vec<u64> = vec![
            17969798196178076256,
            11665707709819759989,
            10951431181515158020,
            3579292038616548363,
            10118509592602682363,
            17443771735705198462,
            13377572520329546867,
            14508682296746081143,
            4893546612470552209,
            6336673525331960373,
            11561838764281975618,
            11059633984043896779,
            9868821492926223124,
            8542944593892234560,
            16032511601890581022,
            7695618831687895843,
            3022120853498371991,
            10947144254984955663,
            9610393059743614339,
            15789984702209978323,
            11507290808992874732,
            8544810361379315625,
            2514505793494337915,
            850311302373851643,
            7471907821583805965,
            4210030402637213022,
            8972286049238152528,
            6681114115741410788,
            2685098934871240750,
            2139759768222668590,
            2799225570783756712,
            11396307301704657495,
            2332683098440312172,
            16891807780001075019,
            16101050408588029653,
            2166370321558883233,
            14347854182549458900,
            14870004927507601360,
            3829254823678205579,
            17782543709877246435,
            7277929836184840438,
            13863079610485812762,
            14662200988771904589,
            6672398562077085777,
            66777234416189871,
            15924117063149907257,
            6229148718437590055,
            1638093424585473056,
            14367533103209980792,
            16637177478581157851,
            12320578874797544608,
            5525203673130288913,
            12739208142136693949,
            7274347089801542680,
            12351602613277061741,
            17125747956453506618,
            10635300702845432897,
            11283900298519168852,
            6267730125728532089,
            6296448508630265012,
            7153662322503197733,
            3465641969518729637,
            5709773413573155800,
            12003525292412621280,
            16270161387413352342,
            940965712129208822,
            5478647835886937857,
            16335137693960490690,
            8252840661133315950,
            6428248950382308744,
            13786474370928231328,
            6201154714147184575,
            16782860777513224723,
            280149148633323228,
            9035407342227835266,
            7068991868375671706,
            2917604801394598957,
            18095856358488458570,
            16227653644574957878,
            17166677276833716419,
            3944277862029338604,
            612646692548617085,
            4310549958317174570,
            12953000510153663529,
            17300537809948341876,
            17579915490199390939,
            14003861546892046105,
            5739963583382643213,
            9985440481510871786,
            9248904368963792352,
            16372336062647489582,
            3041120517596926580,
            534169420124732846,
            4380305661707124107,
            6086724823717347500,
            4223812652099977742,
            11821940997188852439,
            3113343199335373250,
            12170239841078256671,
            17601991461376493256,
            2719634618982556897,
            13195032244809056548,
            10940064296383414607,
            13118802247197706107,
            9565873339643899713,
            14498228386269699899,
            14303271799080240203,
            10557583455300725168,
            5379265057433961705,
            9910550943040392923,
            17035796941001874549,
            7588103056536899678,
            6854744857749589798,
            16544601384842310485,
            3425357872641039597,
            16942896329766599562,
            8845795317407401173,
            17878969833496333473,
            13171171692927769339,
            14950090815969900836,
            12964563090615706228,
            4703872336418279115,
            13577359609057966785,
            5192177735187679108,
            6292649002597521144,
            6149438598104247946,
            14280885284356516174,
            2712424534720961495,
        ];

        let fixed_rng = FixedRandom::new(rands);
        // let fixed_rng = FixedRandom::new(random_numbers_prob);
        // let fixed_rng = FixedRandom1::new(10203421957101150929);

        let mut sampler =
            TernarySampler::new(fixed_rng, ring, Ternary::Probability(1.0 / 3.0), true).unwrap();
        let pol = sampler.read_new();

        // 这里的 expected_coeffs 需要从 Go 测试输出中获取
        let expected_coeffs: Vec<u64> = vec![
            // 在这里填入 Go 测试中 "polynomial coefficients (Probability)" 的输出
        ];

        println!("polynomial coefficients (Probability): {:?}", pol.coeffs[0]);

        assert_eq!(
            pol.coeffs[0], expected_coeffs,
            "Sampled coefficients do not match Go implementation for probability distribution.\nExpected: {:?}\nGot: {:?}",
            expected_coeffs, pol.coeffs[0]
        );
    }

    // #[test]
    // fn test_ternary_sampler_probability() {
    //     let degree = 16 as usize;
    //     let modulus = 97u64;
    //     let p = 0.3;
    //     let montgomery = false;

    //     let ring = Ring::new(degree as u64, vec![modulus]).unwrap();
    //     let prng = ChaCha20Rng::seed_from_u64(12345);

    //     let mut sampler =
    //         TernarySampler::new(prng, ring, Ternary::Probability(p), montgomery).unwrap();
    //     let pol = sampler.read_new();

    //     // test degree of polynomial
    //     assert_eq!(pol.coeffs[0].len(), degree);

    //     // test coefficients are within bounds
    //     for &coeff in &pol.coeffs[0] {
    //         assert!(coeff == 0 || coeff == 1 || coeff == modulus - 1);
    //     }

    //     // test number of non-zero coefficients
    //     let non_zero_count = pol.coeffs[0].iter().filter(|&&x| x != 0).count();
    //     let expected_non_zero = (degree as f64 * (1.0 - p)).round() as usize;
    //     assert!((non_zero_count as i32 - expected_non_zero as i32).abs() <= 3); // 允许一些偏差
    // }

    // #[test]
    // fn test_ternary_sampler_hamming_weight() {
    //     let degree = 16;
    //     let modulus = 97u64;
    //     let hw = 5;
    //     let montgomery = false;

    //     let ring = Ring::new(degree, vec![modulus]).unwrap();
    //     let prng = ChaCha20Rng::seed_from_u64(12345);

    //     let mut sampler =
    //         TernarySampler::new(prng, ring, Ternary::HammingWeight(hw), montgomery).unwrap();
    //     let pol = sampler.read_new();

    //     // test degree of polynomial
    //     assert_eq!(pol.coeffs[0].len(), degree as usize);

    //     // test coefficients are within bounds
    //     for &coeff in &pol.coeffs[0] {
    //         assert!(coeff == 0 || coeff == 1 || coeff == modulus - 1);
    //     }

    //     // test number of non-zero coefficients
    //     let non_zero_count = pol.coeffs[0].iter().filter(|&&x| x != 0).count();
    //     assert_eq!(non_zero_count, hw);
    // }

    // #[test]
    // fn test_ternary_sampler_consistency_hamming_weight() {
    //     let log_n = 10;
    //     let qi: Vec<u64> = vec![
    //         2305843009137934337,
    //         2305843009132953601,
    //         2305843009131642881,
    //         2305843009127448577,
    //         2305843009117224961,
    //         2305843009114341377,
    //         2305843009110409217,
    //         2305843009097564161,
    //         2305843009090748417,
    //         2305843009090486273,
    //         2305843009081311233,
    //         2305843009071087617,
    //         2305843009069514753,
    //         2305843009063223297,
    //     ];

    //     let ring = Ring::new(1 << log_n, qi).unwrap();

    //     let random_numbers_hw: Vec<u64> = vec![12633818882151087297];

    //     let fixed_rng = FixedRandom::new(random_numbers_hw);
    //     let mut sampler =
    //         TernarySampler::new(fixed_rng, ring, Ternary::HammingWeight(128), false).unwrap();

    //     let pol = sampler.read_new();

    //     // 这里的 expected_coeffs 需要从 Go 测试输出中获取
    //     let expected_coeffs: Vec<u64> = vec![
    //         // 在这里填入 Go 测试中 "Sampled polynomial coefficients (Hamming Weight)" 的输出
    //     ];

    //     assert_eq!(
    //         pol.coeffs[0], expected_coeffs,
    //         "Sampled coefficients do not match Go implementation for Hamming weight distribution.\nExpected: {:?}\nGot: {:?}",
    //         expected_coeffs, pol.coeffs[0]
    //     );

    //     // 验证 Hamming weight
    //     let non_zero_count = pol.coeffs[0].iter().filter(|&&x| x != 0).count();
    //     assert_eq!(
    //         non_zero_count, 128,
    //         "Expected Hamming weight 128, but got {}",
    //         non_zero_count
    //     );
    // }
}

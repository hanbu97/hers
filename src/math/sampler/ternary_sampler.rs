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

            println!(
                "randomBytesCoeffs: first 10 bytes: {:?}",
                &random_bytes_coeffs[..10]
            );
            println!(
                "randomBytesCoeffs: last 10 bytes: {:?}",
                &random_bytes_coeffs[random_bytes_coeffs.len() - 10..]
            );
            println!(
                "randomBytesCoeffs: total length: {}",
                random_bytes_coeffs.len()
            );
            println!(
                "randomBytesSign: first 10 bytes: {:?}",
                &random_bytes_sign[..10]
            );
            println!(
                "randomBytesSign: last 10 bytes: {:?}",
                &random_bytes_sign[random_bytes_sign.len() - 10..]
            );
            println!("randomBytesSign: total length: {}", random_bytes_sign.len());

            for i in 0..n {
                let coeff = (random_bytes_coeffs[i >> 3] >> (i & 7)) & 1;
                let sign = (random_bytes_sign[i >> 3] >> (i & 7)) & 1;
                let index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1);

                for (j, &qi) in moduli.iter().enumerate() {
                    let old_value = pol.coeffs[j][i];
                    let lut_value = self.matrix_values[j][index as usize];
                    let new_value = f(old_value, lut_value, qi);
                    pol.coeffs[j][i] = new_value;

                    if i < 5 {
                        println!("Coefficient {}, Modulus {}:", i, j);
                        println!("  coeff: {}, sign: {}, index: {}", coeff, sign, index);
                        println!(
                            "  oldValue: {}, lutValue: {}, qi: {}",
                            old_value, lut_value, qi
                        );
                        println!("  newValue: {}\n", new_value);
                    }
                }
            }
            // } else {
            //     let mut random_bytes = vec![0u8; n];
            //     self.prng.fill_bytes(&mut random_bytes);

            //     // println!("random_bytes: first 10 bytes: {:?}", &random_bytes[..10]);
            //     // println!("random_bytes: last 10 bytes: {:?}", &random_bytes[n - 10..]);
            //     // println!("random_bytes: total length: {}", random_bytes.len());
            //     println!("random_bytes: {:?}", &random_bytes);

            //     let mut pointer = 0u8;
            //     let mut byte_pointer = 0;

            //     for i in 0..n {
            //         let (coeff, sign, new_bytes, new_pointer, new_byte_pointer) =
            //             self.kysampling(&mut random_bytes, pointer, byte_pointer, n);
            //         random_bytes = new_bytes;
            //         pointer = new_pointer;
            //         byte_pointer = new_byte_pointer;

            //         let index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1);

            //         for (j, &qi) in moduli.iter().enumerate() {
            //             let old_value = pol.coeffs[j][i];
            //             let lut_value = self.matrix_values[j][index as usize];
            //             let new_value = f(old_value, lut_value, qi);
            //             pol.coeffs[j][i] = new_value;

            //             if i < 5 {
            //                 println!("Coefficient {}, Modulus {}:", i, j);
            //                 println!("  coeff: {}, sign: {}, index: {}", coeff, sign, index);
            //                 println!(
            //                     "  oldValue: {}, lutValue: {}, qi: {}",
            //                     old_value, lut_value, qi
            //                 );
            //                 println!("  newValue: {}\n", new_value);
            //             }
            //         }
            //     }
            // }
        } else {
            let mut random_bytes = vec![0u8; n];
            self.prng.fill_bytes(&mut random_bytes);

            let mut pointer = 0u8;
            let mut byte_pointer = 0;

            for i in 0..n {
                let (coeff, sign, new_bytes, new_pointer, new_byte_pointer) =
                    self.kysampling(&mut random_bytes, pointer, byte_pointer, n);
                random_bytes = new_bytes;
                pointer = new_pointer;
                byte_pointer = new_byte_pointer;

                let index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1);

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

    // fn kysampling(
    //     &mut self,
    //     random_bytes: &mut Vec<u8>,
    //     pointer: u8,
    //     byte_pointer: usize,
    //     byte_length: usize,
    // ) -> (u64, u64, Vec<u8>, u8, usize) {
    //     // println!(
    //     //     "Rust kysampling input: randomBytes={:?}, pointer={}, bytePointer={}, byteLength={}",
    //     //     &random_bytes[..10],
    //     //     pointer,
    //     //     byte_pointer,
    //     //     byte_length
    //     // );

    //     let mut d = 0;
    //     let mut col = 0;
    //     let col_len = self.matrix_proba.len();
    //     let mut current_pointer = pointer;
    //     let mut current_byte_pointer = byte_pointer;

    //     loop {
    //         for i in current_pointer..8 {
    //             d = (d << 1) + 1 - ((random_bytes[current_byte_pointer] >> i) & 1) as i32;
    //             // println!("Rust kysampling iteration: i={}, d={}", i, d);

    //             if d > col_len as i32 - 1 {
    //                 println!("Rust kysampling: Recursion due to out of bounds");
    //                 return self.kysampling(random_bytes, i, current_byte_pointer, byte_length);
    //             }

    //             for row in (0..col_len).rev() {
    //                 d -= self.matrix_proba[row][col] as i32;

    //                 if d == -1 {
    //                     let sign = if i == 7 {
    //                         current_pointer = 0;
    //                         current_byte_pointer += 1;
    //                         if current_byte_pointer >= byte_length {
    //                             current_byte_pointer = 0;
    //                             self.prng.fill_bytes(random_bytes);
    //                         }
    //                         random_bytes[current_byte_pointer] & 1
    //                     } else {
    //                         current_pointer = i;
    //                         (random_bytes[current_byte_pointer] >> (i + 1)) & 1
    //                     };

    //                     // println!(
    //                     //     "Rust kysampling output: row={}, sign={}, pointer={}, bytePointer={}",
    //                     //     row,
    //                     //     sign,
    //                     //     current_pointer + 1,
    //                     //     current_byte_pointer
    //                     // );
    //                     return (
    //                         row as u64,
    //                         sign as u64,
    //                         random_bytes.to_vec(),
    //                         current_pointer + 1,
    //                         current_byte_pointer,
    //                     );
    //                 }
    //             }

    //             col += 1;
    //         }

    //         current_pointer = 0;
    //         current_byte_pointer += 1;

    //         if current_byte_pointer >= byte_length {
    //             current_byte_pointer = 0;
    //             self.prng.fill_bytes(random_bytes);
    //         }
    //     }
    // }

    fn kysampling(
        &mut self,
        random_bytes: &mut Vec<u8>,
        pointer: u8,
        byte_pointer: usize,
        byte_length: usize,
    ) -> (u64, u64, Vec<u8>, u8, usize) {
        // println!(
        //     "Entering kysampling: pointer={}, byte_pointer={}",
        //     pointer, byte_pointer
        // );
        let mut d = 0;
        let mut col = 0;
        let col_len = self.matrix_proba.len();
        let mut current_pointer = pointer;
        let mut current_byte_pointer = byte_pointer;

        loop {
            for i in current_pointer..8 {
                d = (d << 1) + 1 - ((random_bytes[current_byte_pointer] >> i) & 1) as i32;
                // println!("d={}, col={}, i={}", d, col, i);

                if d > col_len as i32 - 1 {
                    println!("Resampling due to d > col_len - 1");
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
                                println!("Resampling due to end of buffer");
                                current_byte_pointer = 0;
                                self.prng.fill_bytes(random_bytes);
                            }
                            random_bytes[current_byte_pointer] & 1
                        } else {
                            current_pointer = i;
                            (random_bytes[current_byte_pointer] >> (i + 1)) & 1
                        };

                        // println!("Exiting kysampling: row={}, sign={}", row, sign);
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

            println!("Moved to next byte");
            current_pointer = 0;
            current_byte_pointer += 1;

            if current_byte_pointer >= byte_length {
                println!("Resampling due to end of buffer");
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
            7685662699597550763,
            14953292364217102864,
            7222963262361140478,
            16893673922650194115,
            3793088677188345345,
            10087750606732716664,
            8188246280220118973,
            13027583278162930042,
            3924303183778336658,
            7781126008056482024,
            9266264004449922723,
            4391316801689136344,
            1248329669817329268,
            14223357949599275805,
            16730467406247098097,
            55835883205830326,
            14432165027754048161,
            9378575656944622186,
            13357140779710944788,
            8569939460864962352,
            17155089719481595638,
            8777988096342968242,
            14198734893118370825,
            7376626351645426620,
            13112843769534070146,
            3637401487818469167,
            9729926306050373770,
            11933025923153025783,
            12698395781976798179,
            15123698131689175580,
            18409852985484023396,
            1837766387619311271,
            6827264424594972246,
            7182773218724919137,
            1048705044999262477,
            13094520649124109253,
            8381398349756917253,
            14592029365262335357,
            16191177145908300931,
            1169734034294645125,
            6180306290194958031,
            3977595213254218215,
            4451462615574484456,
            10737617763687039816,
            12216633652915601373,
            11158199973728979367,
            13380240135394664655,
            8637589460057531846,
            2565976213863691446,
            8984110935766527066,
            8534377656116200451,
            15469061575852795875,
            4623261720918801142,
            7781632600742237377,
            16317293225425339208,
            4930547598891212082,
            13184128844119340211,
            10574757875384770907,
            7529893994980538324,
            15065039568021523743,
            4751135521472672086,
            14709279943108410984,
            15805275110074818584,
            7841833307574723207,
            1768914828312838226,
            504393574054717531,
            13577523483749371845,
            18396287344711018237,
            12152887487939839129,
            5822194446982284526,
            16081078526038460591,
            10531111026359760537,
            17304988024404459095,
            7980101475738877500,
            4541906559382745604,
            5340399753148586098,
            12559850761087399508,
            275263686075911941,
            8069126251616463418,
            3962033310158913363,
            6160306820045871335,
            3368659379820411299,
            17624727451012183376,
            11034718923000664038,
            16351446891046417047,
            3317224644097274099,
            5439418419354258000,
            768326203528793613,
            13177569846457396345,
            12451828987082913157,
            12586816336386861518,
            14544402751477173997,
            7407904227387548144,
            2350907731690842234,
            11566588287762536628,
            14798968851113923752,
            8135423025063555165,
            9055218322689678584,
            14645145299249945917,
            6822074724995915919,
            3867712485101722095,
            17014973345428190323,
            5106968525718978643,
            3988764555662818515,
            10166198558002146856,
            8091520563969929568,
            8015263996621946405,
            8254969966022434514,
            7431691329007243351,
            17976893343225990381,
            6816875446023252941,
            12748270278000128234,
            6312797555732080511,
            16972592473627671507,
            2373563836649872390,
            16437549214401327645,
            15466197327020354790,
            5724593706030481080,
            7950003937924764077,
            1006372673938639914,
            1068526744112686129,
            10110141121242311288,
            225826799023088687,
            4991598817260496872,
            7299036963762823581,
            1761750928296649414,
            3833830163377678329,
            3936163432609988722,
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

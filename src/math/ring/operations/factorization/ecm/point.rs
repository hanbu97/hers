use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Debug, Clone)]
pub struct CurveParams {
    pub a_24: BigUint,
    pub modulus: BigUint,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    pub x_cord: BigUint,
    pub z_cord: BigUint,
}

impl Point {
    pub fn new(x_cord: BigUint, z_cord: BigUint) -> Point {
        Point { x_cord, z_cord }
    }

    pub fn add(&self, q: &Point, diff: &Point, params: &CurveParams) -> Point {
        let u = (&self.x_cord + &params.modulus - &self.z_cord) * (&q.x_cord + &q.z_cord)
            % &params.modulus;
        let v = (&self.x_cord + &self.z_cord) * (&q.x_cord + &params.modulus - &q.z_cord)
            % &params.modulus;
        let add = (&u + &v) % &params.modulus;
        let subt = (&u + &params.modulus - &v) % &params.modulus;
        let x_cord = (&diff.z_cord * &add * &add) % &params.modulus;
        let z_cord = (&diff.x_cord * &subt * &subt) % &params.modulus;

        Point::new(x_cord, z_cord)
    }

    pub fn double(&self, params: &CurveParams) -> Point {
        let u = (&self.x_cord + &self.z_cord).modpow(&BigUint::from(2u32), &params.modulus);
        let v = (&self.x_cord + &params.modulus - &self.z_cord)
            .modpow(&BigUint::from(2u32), &params.modulus);
        let diff = (&u + &params.modulus - &v) % &params.modulus;
        let x_cord = (&u * &v) % &params.modulus;
        let z_cord = (&v + &params.a_24 * &diff) * &diff % &params.modulus;

        Point::new(x_cord, z_cord)
    }

    pub fn mont_ladder(&self, k: &BigUint, params: &CurveParams) -> Point {
        let mut q = self.clone();
        let mut r = self.double(params);

        for i in format!("{:b}", k)[1..].chars() {
            if i == '1' {
                q = r.add(&q, self, params);
                r = r.double(params);
            } else {
                r = q.add(&r, self, params);
                q = q.double(params);
            }
        }
        q
    }

    pub fn eq(&self, other: &Self, params: &CurveParams) -> bool {
        let self_inv = Self::mod_inverse(&self.z_cord, &params.modulus).unwrap();
        let other_inv = Self::mod_inverse(&other.z_cord, &params.modulus).unwrap();
        (&self_inv * &self.x_cord) % &params.modulus
            == (&other_inv * &other.x_cord) % &params.modulus
    }

    pub fn mod_inverse(a: &BigUint, m: &BigUint) -> Option<BigUint> {
        let (mut t, mut newt) = (BigUint::zero(), BigUint::one());
        let (mut r, mut newr) = (m.clone(), a.clone());

        while !newr.is_zero() {
            let quotient = &r / &newr;
            (t, newt) = (newt.clone(), t + m - (&quotient * &newt) % m);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }

        if r > BigUint::one() {
            None
        } else {
            Some(t % m)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;

    fn create_test_params() -> CurveParams {
        CurveParams {
            a_24: 7.to_biguint().unwrap(),
            modulus: 29.to_biguint().unwrap(),
        }
    }

    #[test]
    fn test_point_add() {
        let params = create_test_params();
        let p1 = Point::new(11.to_biguint().unwrap(), 16.to_biguint().unwrap());
        let p2 = Point::new(13.to_biguint().unwrap(), 10.to_biguint().unwrap());
        let p3 = p2.add(&p1, &p1, &params);

        assert_eq!(p3.x_cord, 23.to_biguint().unwrap());
        assert_eq!(p3.z_cord, 17.to_biguint().unwrap());
    }

    #[test]
    fn test_point_double() {
        let params = create_test_params();
        let p1 = Point::new(11.to_biguint().unwrap(), 16.to_biguint().unwrap());
        let p2 = p1.double(&params);

        assert_eq!(p2.x_cord, 13.to_biguint().unwrap());
        assert_eq!(p2.z_cord, 10.to_biguint().unwrap());
    }

    #[test]
    fn test_point_mont_ladder() {
        let params = create_test_params();
        let p1 = Point::new(11.to_biguint().unwrap(), 16.to_biguint().unwrap());
        let p3 = p1.mont_ladder(&3.to_biguint().unwrap(), &params);

        assert_eq!(p3.x_cord, 23.to_biguint().unwrap());
        assert_eq!(p3.z_cord, 17.to_biguint().unwrap());
    }

    #[test]
    fn test_point() {
        let modulus = 101.to_biguint().unwrap();
        let a: BigUint = 10.to_biguint().unwrap();
        let a_24: BigUint = (&a + 2.to_biguint().unwrap())
            * Point::mod_inverse(&4.to_biguint().unwrap(), &modulus).unwrap()
            % &modulus;

        let params = CurveParams {
            a_24,
            modulus: modulus.clone(),
        };

        let p1 = Point::new(10.to_biguint().unwrap(), 17.to_biguint().unwrap());
        let p2 = p1.double(&params);
        assert!(p2.eq(
            &Point::new(68.to_biguint().unwrap(), 56.to_biguint().unwrap()),
            &params
        ));

        let p4 = p2.double(&params);
        assert!(p4.eq(
            &Point::new(22.to_biguint().unwrap(), 64.to_biguint().unwrap()),
            &params
        ));

        let p8 = p4.double(&params);
        assert!(p8.eq(
            &Point::new(71.to_biguint().unwrap(), 95.to_biguint().unwrap()),
            &params
        ));

        let p16 = p8.double(&params);
        assert!(p16.eq(
            &Point::new(5.to_biguint().unwrap(), 16.to_biguint().unwrap()),
            &params
        ));

        let p32 = p16.double(&params);
        assert!(p32.eq(
            &Point::new(33.to_biguint().unwrap(), 96.to_biguint().unwrap()),
            &params
        ));

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1, &params);
        assert!(p3.eq(
            &Point::new(1.to_biguint().unwrap(), 61.to_biguint().unwrap()),
            &params
        ));

        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1, &params);
        assert!(p5.eq(
            &Point::new(49.to_biguint().unwrap(), 90.to_biguint().unwrap()),
            &params
        ));
        assert!(p5.eq(&p4.add(&p1, &p3, &params), &params));

        // p6 = 2*p3
        let p6 = p3.double(&params);
        assert!(p6.eq(
            &Point::new(87.to_biguint().unwrap(), 43.to_biguint().unwrap()),
            &params
        ));
        assert!(p6.eq(&p4.add(&p2, &p2, &params), &params));

        // p7 = p5 + p2
        let p7 = p5.add(&p2, &p3, &params);
        assert!(p7.eq(
            &Point::new(69.to_biguint().unwrap(), 23.to_biguint().unwrap()),
            &params
        ));
        assert!(p7.eq(&p4.add(&p3, &p1, &params), &params));
        assert!(p7.eq(&p6.add(&p1, &p5, &params), &params));

        // p9 = p5 + p4
        let p9 = p5.add(&p4, &p1, &params);
        assert!(p9.eq(
            &Point::new(56.to_biguint().unwrap(), 99.to_biguint().unwrap()),
            &params
        ));
        assert!(p9.eq(&p6.add(&p3, &p3, &params), &params));
        assert!(p9.eq(&p7.add(&p2, &p5, &params), &params));
        assert!(p9.eq(&p8.add(&p1, &p7, &params), &params));

        assert!(p5.eq(&p1.mont_ladder(&5.to_biguint().unwrap(), &params), &params));
        assert!(p9.eq(&p1.mont_ladder(&9.to_biguint().unwrap(), &params), &params));
        assert!(p16.eq(&p1.mont_ladder(&16.to_biguint().unwrap(), &params), &params));
        assert!(p9.eq(&p3.mont_ladder(&3.to_biguint().unwrap(), &params), &params));
    }
}

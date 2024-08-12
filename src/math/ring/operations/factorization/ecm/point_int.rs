use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

#[derive(Debug, Clone, Default)]
pub struct Point {
    pub x_cord: BigInt,
    pub z_cord: BigInt,
    pub a_24: BigInt,
    pub modulus: BigInt,
}

impl Point {
    pub fn new(x_cord: BigInt, z_cord: BigInt, a_24: BigInt, modulus: BigInt) -> Point {
        Point {
            x_cord,
            z_cord,
            a_24,
            modulus,
        }
    }

    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let u = (&self.x_cord - &self.z_cord) * (&q.x_cord + &q.z_cord);
        let v = (&self.x_cord + &self.z_cord) * (&q.x_cord - &q.z_cord);
        let add = &u + &v;
        let subt = &u - &v;
        let x_cord = (&diff.z_cord * &add * &add).mod_floor(&self.modulus);
        let z_cord = (&diff.x_cord * &subt * &subt).mod_floor(&self.modulus);

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    pub fn double(&self) -> Point {
        let u = (&self.x_cord + &self.z_cord).pow(2u32);
        let v = (&self.x_cord - &self.z_cord).pow(2u32);
        let diff = &u - &v;
        let x_cord = (&u * &v).mod_floor(&self.modulus);
        let z_cord = ((&v + &self.a_24 * &diff) * &diff).mod_floor(&self.modulus);

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    pub fn mont_ladder(&self, k: &BigInt) -> Point {
        let mut q = self.clone();
        let mut r = self.double();

        for i in format!("{:b}", k)[1..].chars() {
            if i == '1' {
                q = r.add(&q, self);
                r = r.double();
            } else {
                r = q.add(&r, self);
                q = q.double();
            }
        }
        q
    }

    fn mod_inverse(a: &BigInt, m: &BigInt) -> Option<BigInt> {
        let (mut t, mut newt) = (BigInt::zero(), BigInt::one());
        let (mut r, mut newr) = (m.clone(), a.clone());

        while !newr.is_zero() {
            let quotient = &r / &newr;
            (t, newt) = (newt.clone(), t - &quotient * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }

        if r > BigInt::one() {
            None
        } else {
            Some((t.mod_floor(m) + m).mod_floor(m))
        }
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        if self.a_24 != other.a_24 || self.modulus != other.modulus {
            false
        } else {
            let self_inv = Point::mod_inverse(&self.z_cord, &self.modulus).unwrap();
            let other_inv = Point::mod_inverse(&other.z_cord, &self.modulus).unwrap();
            (&self_inv * &self.x_cord).mod_floor(&self.modulus)
                == (&other_inv * &other.x_cord).mod_floor(&self.modulus)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_point_add() {
        let p1 = Point::new(
            11.to_bigint().unwrap(),
            16.to_bigint().unwrap(),
            7.to_bigint().unwrap(),
            29.to_bigint().unwrap(),
        );
        let p2 = Point::new(
            13.to_bigint().unwrap(),
            10.to_bigint().unwrap(),
            7.to_bigint().unwrap(),
            29.to_bigint().unwrap(),
        );
        let p3 = p2.add(&p1, &p1);

        assert_eq!(p3.x_cord, 23.to_bigint().unwrap());
        assert_eq!(p3.z_cord, 17.to_bigint().unwrap());
    }

    #[test]
    fn test_point_double() {
        let p1 = Point::new(
            11.to_bigint().unwrap(),
            16.to_bigint().unwrap(),
            7.to_bigint().unwrap(),
            29.to_bigint().unwrap(),
        );
        let p2 = p1.double();

        assert_eq!(p2.x_cord, 13.to_bigint().unwrap());
        assert_eq!(p2.z_cord, 10.to_bigint().unwrap());
    }

    #[test]
    fn test_point_mont_ladder() {
        let p1 = Point::new(
            11.to_bigint().unwrap(),
            16.to_bigint().unwrap(),
            7.to_bigint().unwrap(),
            29.to_bigint().unwrap(),
        );
        let p3 = p1.mont_ladder(&3.to_bigint().unwrap());

        assert_eq!(p3.x_cord, 23.to_bigint().unwrap());
        assert_eq!(p3.z_cord, 17.to_bigint().unwrap());
    }

    #[test]
    fn test_point() {
        let modulus = 101.to_bigint().unwrap();
        let a: BigInt = 10.to_bigint().unwrap();
        let a_24: BigInt = (&a + 2.to_bigint().unwrap())
            * Point::mod_inverse(&4.to_bigint().unwrap(), &modulus).unwrap()
            % &modulus;

        let p1 = Point::new(
            10.to_bigint().unwrap(),
            17.to_bigint().unwrap(),
            a_24.clone(),
            modulus.clone(),
        );
        let p2 = p1.double();
        assert_eq!(
            p2,
            Point::new(
                68.to_bigint().unwrap(),
                56.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        let p4 = p2.double();
        assert_eq!(
            p4,
            Point::new(
                22.to_bigint().unwrap(),
                64.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        let p8 = p4.double();
        assert_eq!(
            p8,
            Point::new(
                71.to_bigint().unwrap(),
                95.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        let p16 = p8.double();
        assert_eq!(
            p16,
            Point::new(
                5.to_bigint().unwrap(),
                16.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        let p32 = p16.double();
        assert_eq!(
            p32,
            Point::new(
                33.to_bigint().unwrap(),
                96.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1);
        assert_eq!(
            p3,
            Point::new(
                1.to_bigint().unwrap(),
                61.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1);
        assert_eq!(
            p5,
            Point::new(
                49.to_bigint().unwrap(),
                90.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        assert_eq!(p5, p4.add(&p1, &p3));
        // p6 = 2*p3
        let p6 = p3.double();
        assert_eq!(
            p6,
            Point::new(
                87.to_bigint().unwrap(),
                43.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        assert_eq!(p6, p4.add(&p2, &p2));
        // p7 = p5 + p2
        let p7 = p5.add(&p2, &p3);
        assert_eq!(
            p7,
            Point::new(
                69.to_bigint().unwrap(),
                23.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        assert_eq!(p7, p4.add(&p3, &p1));
        assert_eq!(p7, p6.add(&p1, &p5));
        // p9 = p5 + p4
        let p9 = p5.add(&p4, &p1);
        assert_eq!(
            p9,
            Point::new(
                56.to_bigint().unwrap(),
                99.to_bigint().unwrap(),
                a_24.clone(),
                modulus.clone()
            )
        );
        assert_eq!(p9, p6.add(&p3, &p3));
        assert_eq!(p9, p7.add(&p2, &p5));
        assert_eq!(p9, p8.add(&p1, &p7));

        assert_eq!(p5, p1.mont_ladder(&5.to_bigint().unwrap()));
        assert_eq!(p9, p1.mont_ladder(&9.to_bigint().unwrap()));
        assert_eq!(p16, p1.mont_ladder(&16.to_bigint().unwrap()));
        assert_eq!(p9, p3.mont_ladder(&3.to_bigint().unwrap()));
    }
}

use crate::decode::{decode_f32, decode_f64};
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::f32::consts::LN_2;
use std::ops::{Add, Div, Mul, MulAssign, Neg, Sub};

const ONE: BigExpFloat = BigExpFloat { exp: 0, float: 1.0 };
const ZERO: BigExpFloat = BigExpFloat { exp: 0, float: 0.0 };

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct BigExpFloat {
    exp: i32,
    float: f32,
}

impl BigExpFloat {
    pub fn from_f32(f: f32) -> Self {
        let (zeroed_exp_f, exp) = decode_f32(f);
        BigExpFloat {
            float: zeroed_exp_f,
            exp,
        }
    }

    pub fn from_f64(f: f64) -> Self {
        let (zeroed_exp_f, exp) = decode_f64(f);
        BigExpFloat {
            float: zeroed_exp_f,
            exp,
        }
    }

    pub fn ln(&self) -> Self {
        let (zeroed_exp_f, exp) = decode_f32(self.float.ln() + (self.exp as f32 * LN_2));
        BigExpFloat {
            float: zeroed_exp_f,
            exp,
        }
    }

    pub fn exp(&self) -> Self {
        let base = BigExpFloat::from_f32(self.float.exp());
        if self.exp.is_positive() {
            let mut acc = base;
            for _ in 0..self.exp {
                acc = acc.square()
            }
            acc
        } else if self.exp.is_negative() {
            let mut acc = base;
            for _ in 0..self.exp.neg() {
                acc = acc.sqrt()
            }
            acc
        } else {
            base
        }
    }

    pub fn sqrt(&self) -> Self {
        if self.exp % 2 == 0 {
            let (zeroed_exp_f, exp) = decode_f32(self.float.sqrt());
            let exp = (self.exp / 2) + exp;
            BigExpFloat {
                float: zeroed_exp_f,
                exp,
            }
        } else {
            let (zeroed_exp_f, exp) = decode_f32(self.float.sqrt() * 2.0_f32.sqrt());
            let exp = ((self.exp - 1) / 2) + exp;
            BigExpFloat {
                float: zeroed_exp_f,
                exp,
            }
        }
    }

    pub fn square(&self) -> Self {
        *self * *self
    }

    pub fn as_f64(&self) -> f64 {
        self.float as f64 * 2.0_f64.powi(self.exp)
    }

    pub fn to_tuple(self) -> (i32, f32) {
        (self.exp, self.float)
    }
}

impl Mul for BigExpFloat {
    type Output = BigExpFloat;

    fn mul(self, rhs: Self) -> Self::Output {
        let res = self.float * rhs.float;
        let (zeroed_exp_f, exp) = decode_f32(res);
        BigExpFloat {
            float: zeroed_exp_f,
            exp: self.exp + rhs.exp + exp,
        }
    }
}

impl MulAssign for BigExpFloat {
    fn mul_assign(&mut self, rhs: Self) {
        let res = self.float * rhs.float;
        let (zeroed_exp_f, exp) = decode_f32(res);
        self.float = zeroed_exp_f;
        self.exp += rhs.exp + exp;
    }
}

impl Div for BigExpFloat {
    type Output = BigExpFloat;

    fn div(self, rhs: Self) -> Self::Output {
        let res = self.float / rhs.float;
        let (zeroed_exp_f, exp) = decode_f32(res);
        BigExpFloat {
            float: zeroed_exp_f,
            exp: self.exp - rhs.exp + exp,
        }
    }
}

impl Add for BigExpFloat {
    type Output = BigExpFloat;

    fn add(self, rhs: Self) -> Self::Output {
        if self.exp == rhs.exp {
            let res = self.float + rhs.float;
            let (zeroed_exp_f, exp) = decode_f32(res);
            BigExpFloat {
                float: zeroed_exp_f,
                exp: self.exp + exp,
            }
        } else {
            // Need to normalize the exponents so that the floats are comparable
            let left_minus_right = self.exp - rhs.exp;
            if left_minus_right > 0 {
                let res = self.float + (rhs.float * 2.0_f32.powi(left_minus_right.neg()));
                let (zeroed_exp_f, exp) = decode_f32(res);
                BigExpFloat {
                    float: zeroed_exp_f,
                    exp: self.exp + exp,
                }
            } else {
                // If left_minus_right < 0
                let res = (self.float * 2.0_f32.powi(left_minus_right)) + rhs.float;
                let (zeroed_exp_f, exp) = decode_f32(res);
                BigExpFloat {
                    float: zeroed_exp_f,
                    exp: rhs.exp + exp,
                }
            }
        }
    }
}

impl Sub for BigExpFloat {
    type Output = BigExpFloat;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.exp == rhs.exp {
            let res = self.float - rhs.float;
            let (zeroed_exp_f, exp) = decode_f32(res);
            BigExpFloat {
                float: zeroed_exp_f,
                exp: self.exp + exp,
            }
        } else {
            // Need to normalize the exponents so that the floats are comparable
            let left_minus_right = self.exp - rhs.exp;
            if left_minus_right > 0 {
                let res = self.float - (rhs.float * 2.0_f32.powi(left_minus_right.neg()));
                let (zeroed_exp_f, exp) = decode_f32(res);
                BigExpFloat {
                    float: zeroed_exp_f,
                    exp: self.exp + exp,
                }
            } else {
                // If left_minus_right < 0
                let res = (self.float * 2.0_f32.powi(left_minus_right)) - rhs.float;
                let (zeroed_exp_f, exp) = decode_f32(res);
                BigExpFloat {
                    float: zeroed_exp_f,
                    exp: rhs.exp + exp,
                }
            }
        }
    }
}

impl Neg for BigExpFloat {
    type Output = BigExpFloat;

    fn neg(self) -> Self::Output {
        BigExpFloat {
            exp: self.exp,
            float: self.float.neg(),
        }
    }
}

impl One for BigExpFloat {
    fn one() -> Self {
        ONE
    }
}

impl Zero for BigExpFloat {
    fn zero() -> Self {
        ZERO
    }

    fn is_zero(&self) -> bool {
        if self.exp == 0 && self.float.is_zero() {
            true
        } else {
            false
        }
    }
}

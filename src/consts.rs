use crate::big_exp_float::BigExpFloat;
use serde::{Deserialize, Serialize};

pub const CANONICAL: bool = true;
pub const DEFAULT_K: usize = 15;
pub const DEFAULT_S: usize = 9;
pub const DEFAULT_T: usize = 2;
pub const DEFAULT_TOTAL_KMERS: usize = 63439162;
pub const TOTAL_CANONICAL_15MERS: usize = 536870912;

pub const REF_SUBDIR: &str = "skim";

#[derive(Serialize, Deserialize)]
pub struct BinomialConsts {
    pub(crate) gamma_r: f64,
    pub(crate) gamma_dk: Vec<BigExpFloat>,
    pub(crate) ln_pi: BigExpFloat,
    pub(crate) ln_2_sqrt_e_over_pi: BigExpFloat,
}

impl BinomialConsts {
    pub fn new() -> Self {
        BinomialConsts {
            gamma_r: 10.900511,
            gamma_dk: vec![
                BigExpFloat::from_f64(2.48574089138753565546e-5),
                BigExpFloat::from_f64(1.05142378581721974210),
                BigExpFloat::from_f64(-3.45687097222016235469),
                BigExpFloat::from_f64(4.51227709466894823700),
                BigExpFloat::from_f64(-2.98285225323576655721),
                BigExpFloat::from_f64(1.05639711577126713077),
                BigExpFloat::from_f64(-1.95428773191645869583e-1),
                BigExpFloat::from_f64(1.70970543404441224307e-2),
                BigExpFloat::from_f64(-5.71926117404305781283e-4),
                BigExpFloat::from_f64(4.63399473359905636708e-6),
                BigExpFloat::from_f64(-2.71994908488607703910e-9),
            ],
            ln_pi: BigExpFloat::from_f64(1.1447298858494001741434273513530587116472948129153),
            ln_2_sqrt_e_over_pi: BigExpFloat::from_f64(
                0.6207822376352452223455184457816472122518527279025978,
            ),
        }
    }
}

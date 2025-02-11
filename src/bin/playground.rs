use std::time::Instant;

use itertools::Itertools;
use num_traits::Zero;
use rand::distr::{Distribution, Uniform};
use skim::{
    big_exp_float::BigExpFloat, binomial_sf::sf, consts::BinomialConsts,
    tracing::start_skim_tracing_subscriber,
};
use statrs::distribution::{Binomial, DiscreteCDF};
use tracing::debug;

fn main() {
    start_skim_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    let assembly_count = 50_000;
    let n = 180;

    let dist = Uniform::new(0.0, 0.01).unwrap();
    let mut rng = rand::rng();
    let file_probabilities = dist
        .sample_iter(&mut rng)
        .take(assembly_count)
        .collect_vec();
    let mut pre_calculated = vec![BigExpFloat::zero(); assembly_count * n];
    let consts = BinomialConsts::new();

    let time = Instant::now();
    pre_calculated
        .iter_mut()
        .enumerate()
        .for_each(|(index, orig)| {
            let (file_num, x) = (index / n, (index % n) as u64);
            let n = n as u64;
            let p = file_probabilities[file_num];
            let prob_f64 = Binomial::new(p, n).unwrap().sf(x);

            // If the probability is greater than 0.0, use it
            let prob_big_exp = if prob_f64 > 0.0 {
                BigExpFloat::from_f64(prob_f64)
            } else {
                // Otherwise, compute the probability using big exp
                sf(p, n, x, &consts)
            };

            *orig = prob_big_exp;
        });
    let total_time = time.elapsed().as_secs_f64();
    debug!("total time {} s", total_time);
    debug!(
        "time per computation {}/s",
        (assembly_count as f64 * n as f64) / total_time,
    );

    debug!("first prob value: {}", file_probabilities[0]);
    for (i, f) in pre_calculated[190..200].iter().enumerate() {
        debug!("example {} {:?}", i, f);
    }
}

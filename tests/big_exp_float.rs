use skim::{binomial_sf::sf, consts::BinomialConsts};
use statrs::distribution::{Binomial, DiscreteCDF};

#[test]
fn many_probabilities() {
    let delta = 10.0_f64.powi(-3);

    let consts = BinomialConsts::new();

    let p = 0.05;
    let n = 20;
    let xs = vec![
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    ];

    let binomial = Binomial::new(p, n).unwrap();

    xs.into_iter().for_each(|x| {
        assert!((binomial.sf(x) - sf(p, n, x, &consts).as_f64()).abs() < delta);
    });
}

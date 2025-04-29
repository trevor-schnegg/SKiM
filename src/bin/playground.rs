use skim::{big_exp_float::BigExpFloat, tracing::start_skim_tracing_subscriber};
use tracing::info;

fn main() {
    start_skim_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    // let kmer_len = 15;
    // let syncmer_info = None;

    // info!("start");

    // let total_syncmers = compute_total_kmers(kmer_len, syncmer_info);

    // info!("total syncmers: {}", total_syncmers);

    let test = 10.0_f64.powi(-9);
    let big_exp_float = BigExpFloat::from_f64(test);

    let exp = 50.0_f64;

    info!("big exp float: {:?}", big_exp_float.powf(exp));
    info!("real float: {:e}", test.powf(exp))
}

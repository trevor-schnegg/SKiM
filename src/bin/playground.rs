use skim::{tracing::start_skim_tracing_subscriber, utility::compute_total_kmers};
use tracing::info;

fn main() {
    start_skim_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    let kmer_len = 15;
    let syncmer_info = None;

    info!("start");

    let total_syncmers = compute_total_kmers(kmer_len, syncmer_info);

    info!("total syncmers: {}", total_syncmers);
}

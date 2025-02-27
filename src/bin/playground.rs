use indicatif::ProgressIterator;
use itertools::Itertools;
use skim::tracing::start_skim_tracing_subscriber;
use std::cmp::min;
use tracing::info;

fn reverse_compliment(kmer: usize, kmer_length: usize, kmer_mask: usize) -> usize {
    let mut buffer = 0;
    let mut complement_kmer = (!kmer) & kmer_mask;
    for _ in 0..kmer_length {
        // Pop the right-most letter
        let letter = complement_kmer & 3;
        complement_kmer >>= 2;
        // Add to the right of the buffer
        buffer <<= 2;
        buffer |= letter;
    }
    buffer
}

fn is_syncmer(kmer: usize, kmer_syncmer_diff: usize, smer_mask: usize) -> bool {
    if kmer_syncmer_diff == 0 {
        true
    } else {
        let minimum_index = (0..=kmer_syncmer_diff)
            .map(|i| (kmer >> ((kmer_syncmer_diff - i) << 1)) & smer_mask)
            .position_min()
            .expect("impossible case");

        if minimum_index == 2 {
            true
        } else {
            false
        }
    }
}

fn main() {
    start_skim_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    let kmer_length = 15_usize;
    let syncmer_length = 9_usize;

    let total_kmers = 4_usize.pow(kmer_length as u32);
    let kmer_syncmer_diff = kmer_length - syncmer_length;
    let kmer_mask = 2_usize.pow((kmer_length * 2) as u32) - 1;
    let syncmer_mask = 2_usize.pow((syncmer_length * 2) as u32) - 1;

    info!("start");

    let total_syncmers = (0..total_kmers)
        .progress()
        .filter_map(|kmer| {
            let canonical_kmer = min(kmer, reverse_compliment(kmer, kmer_length, kmer_mask));
            if kmer == canonical_kmer {
                if is_syncmer(kmer, kmer_syncmer_diff, syncmer_mask) {
                    Some(kmer)
                } else {
                    None
                }
            } else {
                None
            }
        })
        .count();

    info!("total syncmers: {}", total_syncmers);
}

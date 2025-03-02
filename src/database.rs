use num_traits::{One, Zero};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::{collections::HashMap, time::Instant, u16, u32};
use tracing::{debug, info};

use crate::{
    big_exp_float::BigExpFloat,
    binomial_sf::sf,
    consts::BinomialConsts,
    kmer_iter::CanonicalKmerIter,
    rle::{
        Block, BlockIter, NaiveRunLengthEncoding, RunLengthEncoding, MAX_RUN, MAX_UNCOMPRESSED_BITS,
    },
    utility::compute_total_kmers,
};

#[derive(Serialize, Deserialize)]
pub struct Database {
    consts: BinomialConsts,
    files: Box<[String]>,
    kmer_len: usize,
    kmer_to_rle_index: HashMap<u32, u32>,
    p_values: Box<[f64]>,
    rles: Box<[RunLengthEncoding]>,
    syncmer_info: Option<(usize, usize)>,
    tax_ids: Box<[usize]>,
}

impl Database {
    pub fn num_files(&self) -> usize {
        self.files.len()
    }

    pub fn from(
        file_bitmaps: Vec<RoaringBitmap>,
        files: Vec<String>,
        tax_ids: Vec<usize>,
        kmer_len: usize,
        syncmer_info: Option<(usize, usize)>,
    ) -> Self {
        let total_kmers = compute_total_kmers(kmer_len, syncmer_info);
        debug!("{} total possible k-mers", total_kmers);

        // Calculate probability of success (p) for each file with a debug logging step in
        // the middle
        let bitmap_sizes = file_bitmaps
            .par_iter()
            .map(|bitmap| bitmap.len())
            .collect::<Vec<u64>>();
        debug!("total bits set: {}", bitmap_sizes.iter().sum::<u64>());
        let p_values = bitmap_sizes
            .into_par_iter()
            .map(|size| size as f64 / total_kmers as f64)
            .collect::<Box<[f64]>>();

        // Initialize the naive RLEs to be the maximum possible size
        let mut kmer_to_naive_rle =
            vec![NaiveRunLengthEncoding::new(); 4_usize.pow(kmer_len as u32)];

        // Construct all naive kmer RLEs from the bitmaps
        info!("constructing naive runs...");
        for (index, bitmap) in file_bitmaps.into_iter().enumerate() {
            for kmer in bitmap {
                kmer_to_naive_rle[kmer as usize].push(index);
            }
        }

        // Filter out naive rles that do not have any blocks
        let filtered_kmers_and_rles = kmer_to_naive_rle
            .into_iter()
            .enumerate()
            .filter_map(|(kmer, naive_rle)| {
                if naive_rle.num_of_blocks() == 0 {
                    None
                } else {
                    Some((kmer as u32, naive_rle))
                }
            })
            .collect::<Vec<(u32, NaiveRunLengthEncoding)>>();

        // Log information about the number of naive runs
        let naive_run_count = filtered_kmers_and_rles
            .par_iter()
            .map(|(_kmer, naive_rle)| naive_rle.num_of_blocks())
            .sum::<usize>();
        debug!("number of naive rle runs: {}", naive_run_count);

        // Compress the database by allowing the use of uncompressed bit sets
        info!("naive runs constructed! allowing uncompressed bit sets...");
        let kmers_and_rles = filtered_kmers_and_rles
            .into_par_iter()
            .map(|(kmer, naive_rle)| (kmer, naive_rle.to_rle()))
            .collect::<Vec<(u32, RunLengthEncoding)>>();

        // Log information about the number of compressed runs
        let compressed_block_num = kmers_and_rles
            .par_iter()
            .map(|(_kmer, rle)| rle.num_of_blocks())
            .sum::<usize>();
        debug!(
            "number of rle runs after allowing uncompressed bit sets: {}",
            compressed_block_num
        );

        // Create a hashmap over the kmers, indicating where each kmer rle is in the vector
        let mut kmer_to_rle_index = HashMap::with_capacity(kmers_and_rles.len());
        let rles = kmers_and_rles
            .into_iter()
            .enumerate()
            .map(|(index, (kmer, rle))| {
                kmer_to_rle_index.insert(kmer, index as u32);
                rle
            })
            .collect::<Box<[RunLengthEncoding]>>();

        Database {
            consts: BinomialConsts::new(),
            files: files.into_boxed_slice(),
            kmer_len,
            kmer_to_rle_index,
            p_values,
            rles,
            syncmer_info,
            tax_ids: tax_ids.into_boxed_slice(),
        }
    }

    pub fn compute_loookup_table(&self, n_max: u64) -> Vec<BigExpFloat> {
        // Including 0 hits, there are n_max + 1 total possible values for the number of hits
        let possible_hit_numbers = (n_max + 1) as usize;

        let mut lookup_table = vec![BigExpFloat::zero(); self.num_files() * possible_hit_numbers];
        lookup_table
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, placeholder_float)| {
                let (file_num, x) = (
                    index / possible_hit_numbers,
                    (index % possible_hit_numbers) as u64,
                );
                let p = self.p_values[file_num];
                let prob_f64 = Binomial::new(p, n_max).unwrap().sf(x);

                // If the probability is greater than 0.0, use it
                let prob_big_exp = if prob_f64 > 0.0 {
                    BigExpFloat::from_f64(prob_f64)
                } else {
                    // Otherwise, compute the probability using big exp
                    sf(p, n_max, x, &self.consts)
                };

                *placeholder_float = prob_big_exp;
            });

        lookup_table
    }

    pub fn lossy_compression(&mut self, compression_level: usize) -> () {
        fn should_compress(compression_level: usize, set_bits: u32, run_reduction: usize) -> bool {
            if run_reduction < 1 {
                false
            } else {
                match compression_level {
                    1 => set_bits == 1 && run_reduction == 2,
                    2 => set_bits == 1 || set_bits == 2 && run_reduction == 2,
                    3 => set_bits == 1 || set_bits == 2 || set_bits <= 4 && run_reduction == 2,
                    _ => false,
                }
            }
        }

        // Since we are frequently interested in unpacking zeros, create a function to do so
        fn get_zeros(run: Option<&Block>) -> Option<usize> {
            match run {
                Some(Block::Zeros(x)) => Some(*x as usize),
                _ => None,
            }
        }
        info!("performing lossy compresseion...");

        let total_set_bits = self
            .rles
            .par_iter()
            .map(|rle| {
                rle.block_iters()
                    .map(|block_iter| match block_iter {
                        BlockIter::Range((start_i, end_i)) => end_i - start_i,
                        BlockIter::BitIter((bit_iter, _start_i)) => bit_iter.count(),
                    })
                    .sum::<usize>()
            })
            .sum::<usize>();
        debug!("total set bits before compression {}", total_set_bits);

        let total_blocks = self
            .rles
            .par_iter()
            .map(|rle| rle.num_of_blocks())
            .sum::<usize>();
        debug!("total blocks before compression {}", total_blocks);

        self.rles.par_iter_mut().for_each(|current_rle| {
            // variable to hold the new lossy compressed blocks as u16s
            let mut compressed_blocks = vec![];

            // peekable iterator over the current runs
            let mut block_iter = current_rle
                .get_raw_blocks()
                .into_iter()
                .map(|block| Block::from_u16(*block))
                .peekable();

            while let Some(curr_block) = block_iter.next() {
                match curr_block {
                    // if the current run is an uncompressed run, lossy compression may apply
                    Block::Uncompressed(bits) => {
                        // the following variables are `Option` types
                        // if `None`, the previous or next run either did not exist or was not a run of zeros
                        // otherwise, it will be `Some(<length of zero run>)`
                        let prev_zeros = get_zeros(compressed_blocks.last());
                        let next_zeros = get_zeros(block_iter.peek());

                        let set_bits = bits.count_ones();

                        // decide what to do based on the previous and next run
                        match (prev_zeros, next_zeros) {
                            (None, None) => {
                                // no lossy compression possible, push and continue
                                compressed_blocks.push(curr_block);
                            }
                            (Some(length), None) => {
                                // try to merge with the previous
                                let blocks_saved =
                                    if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                        1
                                    } else {
                                        0
                                    };
                                if should_compress(compression_level, set_bits, blocks_saved) {
                                    // If we should compress, add the correct number of
                                    // zeros to the run
                                    match compressed_blocks.last_mut().unwrap() {
                                        // intentionally overloading the variable `length` because they are the same thing
                                        Block::Zeros(length) => {
                                            *length += MAX_UNCOMPRESSED_BITS as u16
                                        }
                                        _ => panic!("impossible case"),
                                    }
                                } else {
                                    // if we shouldn't compress, push and continue
                                    compressed_blocks.push(curr_block);
                                }
                            }
                            (None, Some(length)) => {
                                // try to merge with the next
                                let blocks_saved =
                                    if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                        1
                                    } else {
                                        0
                                    };
                                if should_compress(compression_level, set_bits, blocks_saved) {
                                    // get/burn the next run in the iterator which should be zeros
                                    match block_iter.next().unwrap() {
                                        // intentionally overloading the variable `length` because they are the same thing
                                        Block::Zeros(length) => compressed_blocks.push(
                                            Block::Zeros(MAX_UNCOMPRESSED_BITS as u16 + length),
                                        ),
                                        _ => panic!("impossible case"),
                                    }
                                } else {
                                    // if we shouldn't compress, push and continue
                                    compressed_blocks.push(curr_block);
                                }
                            }
                            (Some(prev_length), Some(next_length)) => {
                                // try to merge with both
                                let total_length =
                                    prev_length + MAX_UNCOMPRESSED_BITS + next_length;
                                let blocks_saved = if total_length <= MAX_RUN as usize {
                                    2
                                } else if total_length <= (MAX_RUN as usize) * 2 {
                                    1
                                } else {
                                    0
                                };
                                if should_compress(compression_level, set_bits, blocks_saved) {
                                    match (
                                        compressed_blocks.last_mut().unwrap(),
                                        block_iter.next().unwrap(),
                                    ) {
                                        // once again, intentionally overloading these variables because they are the same
                                        (Block::Zeros(prev_length), Block::Zeros(next_length)) => {
                                            if total_length <= MAX_RUN as usize {
                                                // if the total length can fit as one run, just add it to the existing run
                                                *prev_length +=
                                                    MAX_UNCOMPRESSED_BITS as u16 + next_length;
                                            } else {
                                                // otherwise, max out the previous run and insert the leftover as a new run of zeros
                                                *prev_length = MAX_RUN;
                                                compressed_blocks.push(Block::Zeros(
                                                    (total_length - MAX_RUN as usize) as u16,
                                                ))
                                            }
                                        }
                                        _ => panic!("impossible case"),
                                    }
                                } else {
                                    // if we shouldn't compress, push and continue
                                    compressed_blocks.push(curr_block);
                                }
                            }
                        };
                    }
                    // if the current run isn't an uncompressed run, push it and continue
                    _ => {
                        compressed_blocks.push(curr_block);
                    }
                } // end match
            } // end while

            // At this point, we have created the variable `compressed_blocks` with the lossy compressed blocks in it
            // We can simply assign the new `compressed_blocks` to the original rle
            let raw_compressed_blocks = compressed_blocks
                .into_iter()
                .map(|block| block.to_u16())
                .collect::<Box<[u16]>>();
            *current_rle = RunLengthEncoding::from(raw_compressed_blocks);
        }); // end par_iter_mut/for_each

        // self.rles has now been mutably updated with the requested lossy compression
        let total_set_bits = self
            .rles
            .par_iter()
            .map(|rle| {
                rle.block_iters()
                    .map(|block_iter| match block_iter {
                        BlockIter::Range((start_i, end_i)) => end_i - start_i,
                        BlockIter::BitIter((bit_iter, _start_i)) => bit_iter.count(),
                    })
                    .sum::<usize>()
            })
            .sum::<usize>();
        debug!("total set bits after compression {}", total_set_bits);

        let total_blocks = self
            .rles
            .par_iter()
            .map(|rle| rle.num_of_blocks())
            .sum::<usize>();
        debug!("total blocks after compression {}", total_blocks);

        // Recompute the p_values after
        self.recompute_p_values();
    }

    fn recompute_p_values(&mut self) -> () {
        let total_kmers = compute_total_kmers(self.kmer_len, self.syncmer_info);
        debug!("{} total possible k-mers", total_kmers);

        let mut file2kmer_num = vec![0_usize; self.num_files()];

        self.rles.iter().for_each(|rle| {
            rle.block_iters().for_each(|block_iter| match block_iter {
                BlockIter::BitIter((bit_iter, start_i)) => {
                    bit_iter.map(|i| i + start_i).for_each(|i| {
                        file2kmer_num[i] += 1;
                    });
                }
                BlockIter::Range((start_i, end_i)) => {
                    file2kmer_num[start_i..end_i].iter_mut().for_each(|count| {
                        *count += 1;
                    });
                }
            });
        });

        let p_values = file2kmer_num
            .into_par_iter()
            .map(|kmer_num| kmer_num as f64 / total_kmers as f64)
            .collect::<Box<[f64]>>();

        self.p_values = p_values;
    }

    pub fn classify(
        &self,
        read: &[u8],
        cutoff_threshold: BigExpFloat,
        n_max: u64,
        lookup_table: &Vec<BigExpFloat>,
    ) -> (Option<(&str, usize)>, (f64, f64)) {
        // Create a vector to store the hits
        let mut num_hits = vec![0_u64; self.num_files()];

        // Create a variable to track the total number of kmers queried
        let mut n_total = 0_u64;

        let hit_lookup_start = Instant::now();
        // For each kmer in the read
        for kmer in
            CanonicalKmerIter::from(read, self.kmer_len, self.syncmer_info).map(|k| k as u32)
        {
            // Lookup the RLE and decompress
            if let Some(rle_index) = self.kmer_to_rle_index.get(&kmer) {
                self.rles[*rle_index as usize].block_iters().for_each(
                    |block_iter| match block_iter {
                        BlockIter::BitIter((bit_iter, start_i)) => {
                            bit_iter.map(|i| i + start_i).for_each(|i| {
                                num_hits[i] += 1;
                            });
                        }
                        BlockIter::Range((start_i, end_i)) => {
                            num_hits[start_i..end_i].iter_mut().for_each(|count| {
                                *count += 1;
                            });
                        }
                    },
                );
            }
            // Increment the total number of queries
            n_total += 1;
        }
        let hit_lookup_time = hit_lookup_start.elapsed().as_secs_f64();

        // Classify the hits
        // Would do this using min_by_key but the Ord trait is difficult to implement for float types
        let prob_calc_start = Instant::now();
        let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
        num_hits
            .iter()
            .zip(self.p_values.iter())
            .enumerate()
            .filter_map(|(index, (n_hits, p))| {
                // This check tries to save runtime in practice
                // Only find the probability if the p-value is going to be < 0.5
                if *n_hits as f64 > (n_total as f64 * p) {
                    // Adjust the number of hits (x) and number of queries (n) based on
                    // the maximum number
                    let x = (*n_hits as f64 * n_max as f64 / n_total as f64).round() as u64;
                    let n = n_max;

                    if n == n_max {
                        // If n is the maximum number, lookup the probability
                        let lookup_position = (index * (n_max + 1) as usize) + x as usize;
                        Some((index, lookup_table[lookup_position]))
                    } else {
                        // Otherwise, perform the computation using f64
                        let prob_f64 = Binomial::new(*p, n).unwrap().sf(x);

                        let prob_big_exp = if prob_f64 > 0.0 {
                            // If the probability is greater than 0.0, use it
                            BigExpFloat::from_f64(prob_f64)
                        } else {
                            // Otherwise, compute the probability using big exp
                            sf(*p, n, x, &self.consts)
                        };

                        Some((index, prob_big_exp))
                    }
                } else {
                    // The p-value will be greater than 0.5 (insignificant)
                    // Don't compute or lookup
                    None
                }
            })
            .for_each(|(index, probability)| {
                // For each index that we computed, compare to find the lowest probability
                // If (for whatever reason) two probabilities are the same, this will use the first one
                if probability < lowest_prob {
                    (lowest_prob_index, lowest_prob) = (index, probability);
                }
            });
        let prob_calc_time = prob_calc_start.elapsed().as_secs_f64();

        if lowest_prob < cutoff_threshold {
            (
                Some((
                    &*self.files[lowest_prob_index],
                    self.tax_ids[lowest_prob_index],
                )),
                (hit_lookup_time, prob_calc_time),
            )
        } else {
            (None, (hit_lookup_time, prob_calc_time))
        }
    }
}

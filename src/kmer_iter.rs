use itertools::Itertools;
use std::cmp::min;
use std::slice::Iter;

const COMPLEMENT: [usize; 4] = [3, 2, 1, 0];

fn base2int(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'a' => Some(0),
        b'C' => Some(1),
        b'c' => Some(1),
        b'G' => Some(2),
        b'g' => Some(2),
        b'T' => Some(3),
        b't' => Some(3),
        _ => None,
    }
}

pub struct KmerIter<'a> {
    canonical: bool,
    char_iter: Iter<'a, u8>,
    curr_kmer: usize,
    curr_rev_comp_kmer: usize,
    first_letter_offset: usize,
    initialized: bool,
    kmer_filter_bits: usize,
    kmer_length: usize,
    kmer_syncmer_diff: usize,
    last_returned: Option<usize>,
    syncmer_filter_bits: usize,
    syncmer_offset: usize,
}

impl<'a> KmerIter<'a> {
    pub fn from(
        sequence: &'a [u8],
        kmer_length: usize,
        canonical: bool,
        syncmer_length: usize,
        syncmer_offset: usize,
    ) -> Self {
        assert!(syncmer_length <= kmer_length);
        assert!(syncmer_offset <= kmer_length - syncmer_length);
        KmerIter {
            canonical,
            char_iter: sequence.iter(),
            curr_kmer: 0,
            curr_rev_comp_kmer: 0,
            first_letter_offset: (kmer_length - 1) * 2,
            initialized: false,
            kmer_filter_bits: 2_usize.pow((kmer_length * 2) as u32) - 1,
            kmer_length,
            kmer_syncmer_diff: kmer_length - syncmer_length,
            last_returned: None,
            syncmer_filter_bits: 2_usize.pow((syncmer_length * 2) as u32) - 1,
            syncmer_offset,
        }
    }

    fn init_next_kmer(&mut self) -> Option<usize> {
        let mut buffer = 0;
        let mut position = 0_usize;
        while position < self.kmer_length {
            match self.char_iter.next() {
                None => {
                    // Reached the end of the sequence without a full k-mer
                    return None;
                }
                Some(char) => {
                    match base2int(*char) {
                        Some(bit_representation) => {
                            buffer <<= 2;
                            buffer |= bit_representation;
                            position += 1;
                        }
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            buffer = 0;
                            position = 0;
                        }
                    }
                }
            }
        }
        // Update the current k-mer and reverse compliment k-mer appropriately
        self.curr_kmer = buffer;
        self.curr_rev_comp_kmer = self.reverse_compliment(buffer);

        // Decide if I need to return this k-mer
        if self.curr_is_syncmer() {
            // Get the correct k-mer
            let return_kmer = if self.canonical {
                min(self.curr_kmer, self.curr_rev_comp_kmer)
            } else {
                self.curr_kmer
            };

            // Check if this k-mer is the same as the one I last returned
            match self.last_returned {
                Some(last_returned_kmer) => {
                    if return_kmer == last_returned_kmer {
                        // If it isn't different, then find the next syncmer that is
                        self.find_next_syncmer()
                    } else {
                        // If it is different, then return it
                        self.last_returned = Some(return_kmer);
                        Some(return_kmer)
                    }
                }
                None => {
                    // If no last k-mer was returned, return this syncmer
                    self.last_returned = Some(return_kmer);
                    Some(return_kmer)
                }
            }
        } else {
            // If the current isn't a syncmer, then find the next one
            self.find_next_syncmer()
        }
    }

    fn find_next_syncmer(&mut self) -> Option<usize> {
        while let Some(char) = self.char_iter.next() {
            match base2int(*char) {
                Some(bit_representation) => {
                    // Update the current k-mer
                    self.curr_kmer <<= 2;
                    self.curr_kmer |= bit_representation;
                    self.curr_kmer &= self.kmer_filter_bits;

                    // Update the current reverse compliment k-mer
                    self.curr_rev_comp_kmer >>= 2;
                    self.curr_rev_comp_kmer |=
                        COMPLEMENT[bit_representation] << self.first_letter_offset;

                    // Decide if I need to return this k-mer
                    if self.curr_is_syncmer() {
                        // Get the correct k-mer
                        let return_kmer = if self.canonical {
                            min(self.curr_kmer, self.curr_rev_comp_kmer)
                        } else {
                            self.curr_kmer
                        };

                        // Check if this k-mer is the same as the one I last returned
                        match self.last_returned {
                            Some(last_returned_kmer) => {
                                if return_kmer == last_returned_kmer {
                                    // If it isn't different, then find the next syncmer that is
                                    continue;
                                } else {
                                    // If it is different, then return it
                                    self.last_returned = Some(return_kmer);
                                    return Some(return_kmer);
                                }
                            }
                            None => {
                                // If no last k-mer was returned, return this syncmer
                                self.last_returned = Some(return_kmer);
                                return Some(return_kmer);
                            }
                        }
                    } else {
                        // If the current isn't a syncmer, then find the next one
                        continue;
                    }
                }
                None => {
                    // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                    return self.init_next_kmer();
                }
            }
        }
        // If we exit the while loop we have no next syncmer
        None
    }

    /// Only call this if I already have an actual k-mer
    fn reverse_compliment(&self, kmer: usize) -> usize {
        let mut buffer = 0;
        let mut complement_kmer = (!kmer) & self.kmer_filter_bits;
        for _ in 0..self.kmer_length {
            // Pop the right-most letter
            let letter = complement_kmer & 3;
            complement_kmer >>= 2;
            // Add to the right of the buffer
            buffer <<= 2;
            buffer |= letter;
        }
        buffer
    }

    fn curr_is_syncmer(&self) -> bool {
        if self.kmer_syncmer_diff == 0 {
            true
        } else {
            let kmer = if self.canonical {
                min(self.curr_kmer, self.curr_rev_comp_kmer)
            } else {
                self.curr_kmer
            };

            let minimum_index = (0..=self.kmer_syncmer_diff)
                .map(|i| (kmer >> ((self.kmer_syncmer_diff - i) << 1)) & self.syncmer_filter_bits)
                .position_min()
                .expect("impossible case");

            if minimum_index == self.syncmer_offset {
                true
            } else {
                false
            }
        }
    }

    pub fn get_curr_kmers(&self) -> (usize, usize) {
        (self.curr_kmer, self.curr_rev_comp_kmer)
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            self.init_next_kmer()
        } else {
            self.find_next_syncmer()
        }
    }
}

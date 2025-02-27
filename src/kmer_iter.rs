use itertools::Itertools;
use std::cmp::min;
use std::slice::Iter;

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

pub struct CanonicalKmerIter<'a> {
    char_iter: Iter<'a, u8>,
    curr_kmer: usize,
    curr_rev_comp_kmer: usize,
    initialized: bool,
    kmer_first_letter_offset: usize,
    kmer_length: usize,
    kmer_mask: usize,
    kmer_smer_diff: usize,
    smer_mask: usize,
    syncmer_offset: usize,
    use_syncmers: bool,
}

impl<'a> CanonicalKmerIter<'a> {
    pub fn from(sequence: &'a [u8], kmer_length: usize, syncmers: Option<(usize, usize)>) -> Self {
        match syncmers {
            Some((smer_length, syncmer_offset)) => {
                assert!(smer_length <= kmer_length);
                assert!(syncmer_offset <= kmer_length - smer_length);
                CanonicalKmerIter {
                    char_iter: sequence.iter(),
                    curr_kmer: usize::MAX,
                    curr_rev_comp_kmer: usize::MAX,
                    initialized: false,
                    kmer_first_letter_offset: (kmer_length - 1) << 1,
                    kmer_length,
                    kmer_mask: (1 << (kmer_length << 1)) - 1,
                    kmer_smer_diff: kmer_length - smer_length,
                    smer_mask: (1 << (smer_length << 1)) - 1,
                    syncmer_offset,
                    use_syncmers: true,
                }
            }
            None => CanonicalKmerIter {
                char_iter: sequence.iter(),
                curr_kmer: usize::MAX,
                curr_rev_comp_kmer: usize::MAX,
                initialized: false,
                kmer_first_letter_offset: (kmer_length - 1) << 1,
                kmer_length,
                kmer_mask: (1 << (kmer_length << 1)) - 1,
                kmer_smer_diff: usize::MAX,
                smer_mask: usize::MAX,
                syncmer_offset: usize::MAX,
                use_syncmers: false,
            },
        }
    }

    fn init_next_kmer(&mut self) -> Option<usize> {
        // Define buffers for the new k-mer
        let mut kmer_buffer = 0;
        let mut rev_comp_kmer_buffer = 0;
        let mut num_kmer_bases = 0_usize;
        while num_kmer_bases < self.kmer_length {
            match self.char_iter.next() {
                None => {
                    // Reached the end of the sequence without a full k-mer
                    return None;
                }
                Some(char) => {
                    match base2int(*char) {
                        Some(c) => {
                            kmer_buffer = (kmer_buffer << 2) | c;
                            rev_comp_kmer_buffer |= (3 - c) << (num_kmer_bases << 1);
                            num_kmer_bases += 1;
                        }
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            // Reset and start over
                            kmer_buffer = 0;
                            rev_comp_kmer_buffer = 0;
                            num_kmer_bases = 0;
                        }
                    }
                }
            }
        }
        // Update the current k-mer and reverse compliment k-mer appropriately
        self.curr_kmer = kmer_buffer;
        self.curr_rev_comp_kmer = rev_comp_kmer_buffer;

        if !self.use_syncmers {
            // If not using syncmers, return the canonical k-mer now
            Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
        } else {
            // If using syncmers, we need to compute if this is a syncmer
            let canonical_kmer = min(self.curr_kmer, self.curr_rev_comp_kmer);
            let min_smer_index = self.kmer_smer_diff
                - (0..=self.kmer_smer_diff)
                    .map(|i| (canonical_kmer >> (i << 1)) & self.smer_mask)
                    .position_min()
                    .expect("impossible case");

            if min_smer_index == self.syncmer_offset {
                // If this is a syncmer, return it
                Some(canonical_kmer)
            } else {
                // Otherwise, look for a syncmer using the other function
                self.find_next_kmer()
            }
        }
    }

    fn find_next_kmer(&mut self) -> Option<usize> {
        while let Some(char) = self.char_iter.next() {
            match base2int(*char) {
                Some(c) => {
                    // Update the current k-mer
                    self.curr_kmer = ((self.curr_kmer << 2) | c) & self.kmer_mask;

                    // Update the current reverse compliment k-mer
                    self.curr_rev_comp_kmer =
                        (self.curr_rev_comp_kmer >> 2) | ((3 - c) << self.kmer_first_letter_offset);

                    if !self.use_syncmers {
                        // If not using syncmers, return the canonical k-mer now
                        return Some(min(self.curr_kmer, self.curr_rev_comp_kmer));
                    } else {
                        // If using syncmers, we need to compute if this is a syncmer
                        let canonical_kmer = min(self.curr_kmer, self.curr_rev_comp_kmer);
                        let min_smer_index = self.kmer_smer_diff
                            - (0..=self.kmer_smer_diff)
                                .map(|i| (canonical_kmer >> (i << 1)) & self.smer_mask)
                                .position_min()
                                .expect("impossible case");

                        if min_smer_index == self.syncmer_offset {
                            // If this is a syncmer, return it
                            return Some(canonical_kmer);
                        } else {
                            // Otherwise, continue the while loop
                            continue;
                        }
                    }
                }
                None => {
                    // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                    return self.init_next_kmer();
                }
            }
        }
        // If we exit the while loop we have no next kmer
        None
    }

    pub fn get_curr_kmers(&self) -> (usize, usize) {
        (self.curr_kmer, self.curr_rev_comp_kmer)
    }
}

impl<'a> Iterator for CanonicalKmerIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            self.init_next_kmer()
        } else {
            self.find_next_kmer()
        }
    }
}

use bio::io::{fasta, fastq};
use indicatif::ProgressIterator;
use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::cmp::min;
use std::fs::File;
use std::fs::{self, DirEntry};
use std::io::{BufReader, ErrorKind};
use std::path::Path;
use std::path::PathBuf;
use tracing::{debug, info, warn};

use crate::kmer_iter::CanonicalKmerIter;

pub const XOR_NUMBER: usize = 188_888_881;

fn is_fasta_file(entry: &DirEntry) -> bool {
    let entry_file_name = entry.file_name().to_str().unwrap().to_string();
    entry_file_name.ends_with(".fna")
        || entry_file_name.ends_with(".fasta")
        || entry_file_name.ends_with(".fa")
}

pub fn get_fasta_files(ref_loc: &Path) -> Vec<PathBuf> {
    let dir_content = fs::read_dir(ref_loc).expect(&*format!(
        "could not read reference directory {:?}",
        ref_loc
    ));

    dir_content
        .par_bridge()
        .into_par_iter()
        .filter_map(|dir_entry| match dir_entry {
            Ok(entry) => {
                if entry.path().is_file() && is_fasta_file(&entry) {
                    Some(entry.path())
                } else {
                    warn!(
                        "reference directory entry {:?} not recognized as a fasta file (did not end with '.fna', '.fasta', or '.fa'), skipping...",
                        entry
                    );
                    None
                }
            }
            Err(e) => {
                warn!("error encountered while reading reference directory: {}. trying to continue...", e);
                None
            }
        })
        .collect::<Vec<PathBuf>>()
}

pub fn get_fasta_iter_of_file(file_path: &Path) -> fasta::Records<BufReader<File>> {
    match fasta::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!(
            "unable to get fasta iter of {:?} because of the following error: {}",
            file_path, error
        ),
    }
}

pub fn get_fastq_iter_of_file(file_path: &Path) -> fastq::Records<BufReader<File>> {
    match fastq::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!(
            "unable to get fastq iter of {:?} because of the following error: {}",
            file_path, error
        ),
    }
}

// Creates a single bitmap containing k-mers from all files, if necessary
pub fn create_bitmap(
    file: PathBuf,
    kmer_len: usize,
    syncmers: Option<(usize, usize)>,
) -> RoaringBitmap {
    let mut bitmap = RoaringBitmap::new();
    let mut record_iter = get_fasta_iter_of_file(&file);
    while let Some(Ok(record)) = record_iter.next() {
        for kmer in CanonicalKmerIter::from(record.seq(), kmer_len, syncmers) {
            bitmap.insert(kmer as u32);
        }
    }
    bitmap
}

pub fn compute_total_kmers(kmer_len: usize, syncmers: Option<(usize, usize)>) -> usize {
    let total_kmers = 4_usize.pow(kmer_len as u32);
    let kmer_mask = (1 << (kmer_len << 1)) - 1;
    info!("computing total possible k-mers...");
    match syncmers {
        Some((smer_len, syncmer_offset)) => {
            let smer_mask = (1 << (smer_len << 1)) - 1;
            let kmer_smer_diff = kmer_len - smer_len;
            (0..total_kmers)
                .progress()
                .filter_map(|kmer| {
                    let canonical_kmer = min(kmer, reverse_compliment(kmer, kmer_len, kmer_mask));
                    if kmer == canonical_kmer {
                        if is_syncmer(kmer, kmer_smer_diff, smer_mask, syncmer_offset) {
                            Some(kmer)
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .count()
        }
        None => (0..total_kmers)
            .progress()
            .filter_map(|kmer| {
                let canonical_kmer = min(kmer, reverse_compliment(kmer, kmer_len, kmer_mask));
                if kmer == canonical_kmer {
                    Some(kmer)
                } else {
                    None
                }
            })
            .count(),
    }
}

fn reverse_compliment(kmer: usize, kmer_len: usize, kmer_mask: usize) -> usize {
    let mut buffer = 0;
    let mut complement_kmer = (!kmer) & kmer_mask;
    for _ in 0..kmer_len {
        // Pop the right-most letter
        let letter = complement_kmer & 3;
        complement_kmer >>= 2;
        // Add to the right of the buffer
        buffer = (buffer << 2) | letter;
    }
    buffer
}

fn is_syncmer(kmer: usize, kmer_smer_diff: usize, smer_mask: usize, syncmer_offset: usize) -> bool {
    if kmer_smer_diff == 0 {
        true
    } else {
        let minimum_index = (0..=kmer_smer_diff)
            .map(|i| (kmer >> ((kmer_smer_diff - i) << 1)) & smer_mask)
            .position_min()
            .expect("impossible case");

        if minimum_index == syncmer_offset {
            true
        } else {
            false
        }
    }
}

pub fn create_ref_subdir(ref_subdir: &PathBuf) {
    info!("creating directory {:?} to store split files", ref_subdir);
    match fs::create_dir(ref_subdir) {
        Ok(()) => debug!("dirctory created successfully!"),
        Err(e) => match e.kind() {
            ErrorKind::AlreadyExists => {
                debug!("directory already exists (likely from a previous execution), continuing...")
            }
            e => panic!("{}", e),
        },
    }
}

pub fn split_record(
    record: fasta::Record,
    max_len: usize,
    overlap_len: usize,
) -> Vec<fasta::Record> {
    let half_overlap_len = overlap_len / 2;

    let num_fragments = record.seq().len().div_ceil(max_len);
    let fragment_len = record.seq().len().div_ceil(num_fragments);

    (0..num_fragments)
        .map(|fragment_index| {
            let start = fragment_index * fragment_len;
            let start_w_overlap = if start == 0 {
                0
            } else {
                start - half_overlap_len
            };
            let end_w_overlap = min(start + fragment_len + half_overlap_len, record.seq().len());

            fasta::Record::with_attrs(
                record.id(),
                Some(&*format!("{}..{}", start_w_overlap, end_w_overlap)),
                &record.seq()[start_w_overlap..end_w_overlap],
            )
        })
        .collect::<Vec<fasta::Record>>()
}

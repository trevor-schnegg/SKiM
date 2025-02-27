use bio::io::{fasta, fastq};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::fs::File;
use std::fs::{self, DirEntry};
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;
use tracing::{error, warn};

use crate::kmer_iter::CanonicalKmerIter;

pub const XOR_NUMBER: usize = 188_888_881;

fn is_fasta_file(entry: &DirEntry) -> bool {
    let entry_file_name = entry.file_name().to_str().unwrap().to_string();
    entry_file_name.ends_with(".fna")
        || entry_file_name.ends_with(".fasta")
        || entry_file_name.ends_with(".fa")
}

pub fn get_fasta_files(reference_loc: &Path) -> Vec<PathBuf> {
    let dir_content = fs::read_dir(reference_loc).expect("could not read reference directory");
    dir_content
        .par_bridge()
        .into_par_iter()
        .filter_map(|dir_entry| match dir_entry {
            Ok(entry) => {
                if entry.path().is_file() && is_fasta_file(&entry) {
                    Some(entry.path())
                } else {
                    warn!(
                        "directory entry {:?} did not end with '.fna', '.fasta', or '.fa', skipping...",
                        entry
                    );
                    None
                }
            }
            Err(e) => {
                error!(
                    "error encountered while reading reference directory {:?}",
                    reference_loc
                );
                error!("{}", e);
                warn!("attempting to continue execution...");
                None
            }
        })
        .collect::<Vec<PathBuf>>()
}

pub fn get_fasta_iter_of_file(file_path: &Path) -> fasta::Records<BufReader<File>> {
    match fasta::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!("{}", error),
    }
}

pub fn get_fastq_iter_of_file(file_path: &Path) -> fastq::Records<BufReader<File>> {
    match fastq::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!("{}", error),
    }
}

// Creates a single bitmap containing k-mers from all files, if necessary
pub fn create_bitmap(
    files: Vec<PathBuf>,
    kmer_len: usize,
    syncmers: Option<(usize, usize)>,
) -> RoaringBitmap {
    let mut bitmap = RoaringBitmap::new();
    for file in files {
        let mut record_iter = get_fasta_iter_of_file(&file);
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_len {
                continue;
            }
            for kmer in CanonicalKmerIter::from(record.seq(), kmer_len, syncmers) {
                bitmap.insert(kmer as u32);
            }
        }
    }
    bitmap
}

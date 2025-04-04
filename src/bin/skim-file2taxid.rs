use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use rayon::prelude::*;
use skim::consts::{DEFAULT_K, REF_SUBDIR};
use skim::io::{create_output_file, load_string2taxid, save_fasta_record_to_file};
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::{
    compute_total_kmers, create_ref_subdir, get_fasta_files, get_fasta_iter_of_file, split_record,
};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use tracing::{info, warn};

const MAX_PROB: f64 = 0.1;

fn get_taxid(accession2taxid: &Option<HashMap<String, usize>>, accession: &str) -> usize {
    match accession2taxid {
        Some(accession2taxid) => *accession2taxid.get(accession).expect(&*format!(
            "accession2taxid was provided but did not contain the key: {}",
            accession
        )),
        None => 0,
    }
}

/// Creates a file2taxid (.f2t) file of the form <fasta-file>\t<taxid> for a given a reference location.
/// Splits files that may be too large for use with the given k-mer size in the process.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, verbatim_doc_comment)]
    /// The accession2taxid/seqid2taxid file.
    /// If not provided, all tax ids will be set to 0.
    /// SKiM will still report the file that each read hits to.
    accession2taxid: Option<String>,

    #[arg(short, long, default_value_t = DEFAULT_K)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the file2taxid (.f2t) file.
    /// If a file is provided, the extention '.skim.f2t' is added.
    /// If a directory is provided, 'skim.f2t' will be the file name.
    output_location: String,

    #[arg(short = 'l', long, default_value_t = 1_000_000)]
    /// The total base pairs of overlap if sequences need to be split
    overlap_length: usize,

    #[arg()]
    /// Directory with FASTA file targets of the reference database
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let kmer_len = args.kmer_length;
    let output_loc_path = Path::new(&args.output_location);
    let overlap_len = args.overlap_length;
    let ref_dir_path = Path::new(&args.reference_directory);

    // Create the output file so it errors if a bad output file is provided before computation
    let mut output_writer = BufWriter::new(create_output_file(output_loc_path, "skim.f2t"));

    info!("k-mer length: {}", kmer_len);
    let total_kmers = compute_total_kmers(kmer_len, None);
    let total_len_allowed = (total_kmers as f64 * MAX_PROB).round() as usize;

    // Get the accession2taxid, if one was provided
    let accession2taxid: Option<HashMap<String, usize>> = match args.accession2taxid {
        None => {
            warn!("no accession2taxid was provided - setting all tax ids to 0");
            warn!("please be sure this is intentional");
            None
        }
        Some(accession2taxid) => {
            let accession2taxid_path = Path::new(&accession2taxid);
            info!("reading accession2taxid at {}", accession2taxid);
            Some(HashMap::from_iter(
                load_string2taxid(accession2taxid_path).into_iter(),
            ))
        }
    };

    let ref_subdir = ref_dir_path.join(REF_SUBDIR);
    let ref_subdir_created = Mutex::new(false);

    info!("searching through files in {}", args.reference_directory);
    let file2taxid = get_fasta_files(ref_dir_path)
        .into_par_iter()
        .progress()
        .map(|file| {
            // Collect records from the fasta file
            let records = get_fasta_iter_of_file(&file)
                .map(|record_result| match record_result {
                    Ok(record) => record,
                    Err(e) => panic!(
                        "could not parse the fasta file {:?} due to the following error {}",
                        file, e,
                    ),
                })
                .collect_vec();

            // Get the tax id of the file based on the first record
            let taxid = get_taxid(
                &accession2taxid,
                records
                    .get(0)
                    .expect(&*format!("no records found in fasta file {:?}", file))
                    .id(),
            );

            // Check if there are too many base pairs
            let total_len_of_records = records.iter().map(|r| r.seq().len()).sum::<usize>();
            if total_len_of_records > total_len_allowed {
                // If there are too many base pairs, create a directory to store the split files
                // This lock ensures that only one thread creates the directory and that all other
                // threads wait while doing so
                {
                    let mut ref_subdir_created = ref_subdir_created.lock().unwrap();
                    if !*ref_subdir_created {
                        info!("a file exceeded the max number of base pairs to reliably be used in classification");
                        info!("skim will try to split such files by record");
                        info!("if an individual record is still too long, the record will be split into fragments with an overlap of {} (as specified by the option -l)", overlap_len);
                        create_ref_subdir(&ref_subdir);
                        *ref_subdir_created = true;
                    }
                }

                // Finally split the records as needed and save them in the new directory
                let file_name = file.file_name().unwrap().to_str().unwrap();
                records.into_iter().enumerate().map(|(record_index, record)| {
                    let output_file_base = ref_subdir.join(file_name);

                    if record.seq().len() > total_len_allowed {
                        // If the record is still too long, split into fragments and save to files
                        split_record(record, total_len_allowed, overlap_len).into_iter().enumerate().map(|(record_subindex, new_record)| {
                            let output_file = output_file_base.with_extension(format!("{}.{}", record_index, record_subindex));
                            save_fasta_record_to_file(new_record, &output_file);
                            (output_file, taxid)
                        }).collect::<Vec<(PathBuf, usize)>>()
                    } else {
                        // If the record is short enough, just save to a new file
                        let output_file = output_file_base.with_extension(record_index.to_string());
                        save_fasta_record_to_file(record, &output_file);
                        vec![(output_file, taxid)]
                    }
                }).flatten().collect::<Vec<(PathBuf, usize)>>()
            } else {
                // If there are not too many base pairs, simply return the file and taxid
                vec![(file, taxid)]
            }
        })
        .flatten()
        .collect::<Vec<(PathBuf, usize)>>();

    info!("writing to output file...");

    for (file_path, taxid) in file2taxid {
        // Write the result to the output file
        let file_name = file_path.file_name().unwrap().to_str().unwrap();

        let relative_file_path = match file_path.parent() {
            Some(parent_dir) => {
                if parent_dir.file_name().unwrap().to_str().unwrap() == REF_SUBDIR {
                    &*format!("{}/{}", REF_SUBDIR, file_name)
                } else {
                    file_name
                }
            }
            None => file_name,
        };
        output_writer
            .write(format!("{}\t{}\n", relative_file_path, taxid).as_bytes())
            .expect("could not write to output file");
    }

    output_writer.flush().unwrap();

    info!("done");
}

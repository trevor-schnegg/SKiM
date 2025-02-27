use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use skim::io::{create_output_file, dump_data_to_file, load_string2taxid};
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::create_bitmap;
use std::path::Path;
use tracing::info;

/// Computes the pairwise distance (.pd) matrix (lower triangle) from the input file2taxid
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 15)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the pairwise distance (.pd) file.
    /// If a file is provided, the extention '.skim.pd' is added.
    /// If a directory is provided, 'skim.pd' will be the file name.
    output_location: String,

    #[arg(short, long, default_value_t = 9)]
    /// Length of syncmer to use in the database
    smer_length: usize,

    #[arg(short = 't', long, default_value_t = 3)]
    /// Offset of syncmer to use in the database
    syncmer_offset: usize,

    #[arg()]
    /// The file2taxid (.f2t) file
    file2taxid: String,

    #[arg()]
    /// Directory with FASTA file targets of the reference database
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let kmer_len = args.kmer_length;
    let output_loc_path = Path::new(&args.output_location);
    let ref_dir_path = Path::new(&args.reference_directory);

    let syncmers = if kmer_len == args.smer_length {
        info!(
            "syncmers disabled: k-mer length ({}) is the same as the syncmer length",
            kmer_len
        );
        None
    } else {
        info!(
            "k-mer length: {}, s-mer length: {}, syncmer offset: {}",
            kmer_len, args.smer_length, args.syncmer_offset
        );
        Some((args.smer_length, args.syncmer_offset))
    };

    // Create the output file so it errors if an incorrect output file is provided before computation
    let output_file = create_output_file(output_loc_path, "skim.pd");

    info!("loading file2taxid at {}", args.file2taxid);
    let file2taxid = load_string2taxid(file2taxid_path);

    info!("creating roaring bitmaps for each group...");
    let bitmaps = file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            // Split the files up if they are grouped
            let file_paths = files
                .split("$")
                .map(|file| ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, syncmers)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("roaring bitmaps created, creating distance matrix...");
    let distances = bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .map(|(index_1, bitmap_1)| {
            bitmaps[..=index_1]
                .iter()
                .enumerate()
                .map(|(index_2, bitmap_2)| {
                    if index_1 == index_2 {
                        0
                    } else {
                        let intersection_size = bitmap_1.intersection_len(bitmap_2);
                        // |A| + |B| - (2 * |A & B|)
                        (bitmap_1.len() + bitmap_2.len() - (2 * intersection_size)) as u32
                    }
                })
                .collect::<Vec<u32>>()
        })
        .collect::<Vec<Vec<u32>>>();

    info!("distance matrix completed! outputting to file...");
    dump_data_to_file(&(distances, file2taxid), output_file)
        .expect("could not output distances to file");

    info!("done!");
}

use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use skim::io::{create_output_file, dump_data_to_file, load_data_from_file, load_string2taxid};
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::create_bitmap;
use std::path::Path;
use tracing::info;

/// Computes the lower triangle of a pairwise distance matrix from the input sequences (or sequence groups)
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 15)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output
    /// If a file, '.skim.pd' is added
    /// If a directory, 'skim.pd' will be the file name
    /// Name means: skim, (p)airwise (d)istances
    output_location: String,

    #[arg(short, long, default_value_t = 9)]
    /// Length of syncmer to use in the database
    smer_length: usize,

    #[arg(short = 't', long, default_value_t = 3)]
    /// Offset of syncmer to use in the database
    syncmer_offset: usize,

    #[arg()]
    /// The original pairwise distances file
    distances: String,

    #[arg()]
    /// The file2taxid map file of the new fasta files
    new_file2taxid: String,

    #[arg()]
    /// Directory with fasta files to extend the distances with
    new_reference_directory: String,

    #[arg()]
    /// Directory with fasta files that the distances were created with
    old_reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_path = Path::new(&args.distances);
    let new_file2taxid_path = Path::new(&args.new_file2taxid);
    let kmer_len = args.kmer_length;
    let new_ref_dir_path = Path::new(&args.new_reference_directory);
    let old_ref_dir_path = Path::new(&args.old_reference_directory);
    let output_loc_path = Path::new(&args.output_location);

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

    // Create the output file
    let output_file = create_output_file(output_loc_path, "skim.pd");

    info!("loading pairwise distances at {}", args.distances);
    let (old_distances, old_file2taxid) =
        load_data_from_file::<(Vec<Vec<u32>>, Vec<(String, usize)>)>(distances_path);
    let old_file2taxid_len = old_file2taxid.len();

    info!("loading new file2taxid at {:?}", new_file2taxid_path);
    let new_file2taxid = load_string2taxid(new_file2taxid_path);

    info!("creating bitmaps for the old file2taxid...");
    let old_bitmaps = old_file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            // Split the files up if they are grouped
            let file_paths = files
                .split("$")
                .map(|file| old_ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, syncmers)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!(
        "{} groups need to be added, creating roaring bitmaps for new file2taxid...",
        new_file2taxid.len()
    );
    let new_bitmaps = new_file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| new_ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, syncmers)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("filling out distance matrix...");
    let all_bitmaps = old_bitmaps
        .into_iter()
        .chain(new_bitmaps.into_iter())
        .collect_vec();

    let new_distances = all_bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .filter_map(|(index_1, bitmap_1)| {
            if index_1 < old_file2taxid_len {
                None
            } else {
                Some(
                    all_bitmaps[..=index_1]
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
                        .collect::<Vec<u32>>(),
                )
            }
        })
        .collect::<Vec<Vec<u32>>>();

    info!("combining and outputting to file...");
    let all_file2taxid = old_file2taxid
        .into_iter()
        .chain(new_file2taxid.into_iter())
        .collect_vec();

    let all_distances = old_distances
        .into_iter()
        .chain(new_distances.into_iter())
        .collect_vec();

    dump_data_to_file(&(all_distances, all_file2taxid), output_file).unwrap();

    info!("done!");
}

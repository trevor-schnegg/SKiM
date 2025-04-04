use clap::Parser;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use skim::consts::{DEFAULT_K, DEFAULT_S, DEFAULT_T};
use skim::database::Database;
use skim::io::{create_output_file, dump_data_to_file, load_string2taxid};
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::create_bitmap;
use std::path::Path;
use tracing::info;

/// Creates a skim database (.db) file from a file2taxid (.f2t) file.
/// For significant database size improvement, the file2taxid should be ordered (.o.f2t).
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = DEFAULT_K)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the database (.db) file.
    /// If a file is provided, the extension '.skim.db' is added.
    /// If a directory is provided, 'skim.db' will be the file name.
    output_location: String,

    #[arg(short, long, default_value_t = DEFAULT_S)]
    /// Length of s-mer to use in the database
    smer_length: usize,

    #[arg(short = 't', long, default_value_t = DEFAULT_T)]
    /// Offset of s-mer to create a syncmer for database
    /// 0 indicates no offset (open syncmers)
    syncmer_offset: usize,

    #[arg()]
    /// The file2taxid (.f2t) file. Preferrably ordered (.o.f2t).
    file2taxid: String,

    #[arg()]
    /// Directory with FASTA files targets of the reference database
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let kmer_len = args.kmer_length;
    let file2taxid_path = Path::new(&args.file2taxid);
    let output_loc_path = Path::new(&args.output_location);
    let ref_dir_path = Path::new(&args.reference_directory);

    // Specifically parse how syncmers should be handled
    let syncmer_info = if kmer_len == args.smer_length {
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

    // Create the output file so it errors if a bad output file is provided before computation
    let output_file = create_output_file(output_loc_path, "skim.db");

    // Load the file2taxid ordering
    info!("loading file2taxid at {}", args.file2taxid);
    let (files, tax_ids): (Vec<String>, Vec<usize>) =
        load_string2taxid(file2taxid_path).into_iter().unzip();

    info!("creating roaring bitmaps for each file...");
    let bitmaps = files
        .par_iter()
        .progress()
        .map(|file| create_bitmap(ref_dir_path.join(file), kmer_len, syncmer_info))
        .collect::<Vec<RoaringBitmap>>();

    info!("constructing database...");
    let database = Database::from(bitmaps, files, tax_ids, kmer_len, syncmer_info);

    info!("dumping to file...");
    dump_data_to_file(&database, output_file).expect("could not serialize database to file");

    info!("done!");
}

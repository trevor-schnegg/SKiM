use clap::Parser;
use skim::database::Database;
use skim::io::{create_output_file, dump_data_to_file, load_data_from_file};
use skim::tracing::start_skim_tracing_subscriber;
use std::path::Path;
use tracing::info;

/// Creates a lossy compressed skim datbase (.cdb) file from a skim database (.db) file.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the database (.cdb) file.
    /// If a file is provided, the extension '.skim.cdb' is added.
    /// If a directory is provided, 'skim.cdb' will be the file name.
    output_location: String,

    #[arg()]
    /// Level of compression: one of [1, 2, 3]
    compression_level: usize,

    #[arg()]
    /// The uncompressed database file
    database: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let database_path = Path::new(&args.database);
    let output_loc_path = Path::new(&args.output_location);

    // Assert that the user provided 1, 2, or 3 as the compression level
    let compression_level = if args.compression_level >= 1 && args.compression_level <= 3 {
        args.compression_level
    } else {
        panic!("compression level was not 1, 2, or 3!")
    };

    // Create the output file so it errors if an incorrect output file is provided before computation
    let output_file = create_output_file(output_loc_path, "skim.cdb");

    info!("loading database at {:?}", database_path);
    let mut database = load_data_from_file::<Database>(database_path);

    info!(
        "compressing database using compression level: {}",
        compression_level
    );
    database.lossy_compression(compression_level);

    info!("dumping to file...");
    dump_data_to_file(&database, output_file).expect("could not output database to file");

    info!("done!");
}

use clap::Parser;
use skim::database::Database;
use skim::io::{create_output_file, dump_data_to_file, load_data_from_file};
use skim::tracing::start_skim_tracing_subscriber;
use std::path::Path;
use tracing::info;

/// Creates a skim database (.db) file from a file2taxid (.f2t) file.
/// For significant database size improvement, the file2taxid should be ordered (.o.f2t).
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the database (.db) file.
    /// If a file is provided, the extension '.skim.db' is added.
    /// If a directory is provided, 'skim.db' will be the file name.
    output_location: String,

    #[arg()]
    /// The database (.db/.cdb) file
    database: String,

    #[arg()]
    /// The file to replace the tax id of
    file: String,

    #[arg()]
    /// The replacement tax id
    taxid: usize,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let database_path = Path::new(&args.database);
    let output_loc_path = Path::new(&args.output_location);

    // Create the output file so it errors if a bad output file is provided before computation
    let output_file = create_output_file(output_loc_path, "fixed.skim.db");

    info!("loading database at {:?}", database_path);
    let mut database = load_data_from_file::<Database>(database_path);

    database.update_taxid(args.file, args.taxid);

    info!("dumping to file...");
    dump_data_to_file(&database, output_file).expect("could not serialize database to file");

    info!("done!");
}

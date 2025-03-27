use clap::Parser;
use skim::{
    io::{create_output_file, load_data_from_file},
    order::{greedy_ordering, ordering_statistics},
    tracing::start_skim_tracing_subscriber,
};
use std::{
    io::{BufWriter, Write},
    path::Path,
};
use tracing::{debug, info};

/// Creates an ordered file2taxid (.o.f2t) file based on a pairwise distance (.pd) matrix.
/// This is done such that the total hamming distance of the ordering is as small as possible.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the ordered file2taxid (.o.f2t) file.
    /// If a file is provided, the extention '.skim.o.f2t' is added.
    /// If a directory is provided, 'skim.o.f2t' will be the file name.
    output_location: String,

    #[arg(short, long, default_value_t = 0)]
    /// Start index of the naive shortest path traversal
    start: usize,

    #[arg()]
    /// The pairwise distances (.pd) file
    distances: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_file = Path::new(&args.distances);
    let output_loc_path = Path::new(&args.output_location);

    // Create the output file so it errors if a bad output file is provided before computation
    let mut output_writer = BufWriter::new(create_output_file(output_loc_path, "skim.o.f2t"));

    info!("loading distances at {}", args.distances);
    let (distances, file2taxid) =
        load_data_from_file::<(Vec<Vec<u32>>, Vec<(String, usize)>)>(distances_file);

    info!("distances loaded! finding ordering...");

    // Perform the greedy solution -- no other options for right now
    let greedy_ordering = greedy_ordering(&distances, args.start);
    let (avg_dist, total_dist) = ordering_statistics(&greedy_ordering, &distances);
    debug!("length of tour: {}", total_dist);
    debug!("average distance between files: {}", avg_dist);

    info!("writing to output file...");

    for index in greedy_ordering {
        let (files_string, taxid) = &file2taxid[index];
        output_writer
            .write(format!("{}\t{}\n", *files_string, *taxid).as_bytes())
            .expect("could not write to output file");
    }

    output_writer.flush().unwrap();

    info!("done!");
}

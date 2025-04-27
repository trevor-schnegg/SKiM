use clap::Parser;
use rayon::prelude::*;
use skim::big_exp_float::BigExpFloat;
use skim::database::Database;
use skim::io::{create_output_file, load_data_from_file};
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::get_fastq_iter_of_file;
use std::io::{BufWriter, Write};
use std::ops::Neg;
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;
use tracing::{debug, info, warn};

/// Classifies the input reads using a skim database (.db/.cdb) file.
/// Output is a readid2file (.r2f) mapping, including the taxid for the file if it was provided during database construction.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 9, verbatim_doc_comment)]
    /// The exponent 'e' used in the equation 10^{-e}.
    /// Any calculated p-value below 10^{-e} will result in a classification.
    exp_cutoff: i32,

    #[arg(short, long, default_value_t = 100, verbatim_doc_comment)]
    /// The fixed number of trials to use in the binomial function.
    n_fixed: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the readid2file (.r2f) file.
    /// If a file is provided, the extension '.skim.r2f' is added.
    /// If a directory is provided, 'skim.r2f' will be the file name.
    output_location: String,

    #[arg()]
    /// The database (.db/.cdb) file
    database: String,

    #[arg()]
    /// FASTQ reads file to query
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let cutoff_threshold = BigExpFloat::from_f64(10.0_f64.powi((args.exp_cutoff).neg()));
    let database_path = Path::new(&args.database);
    let output_loc_path = Path::new(&args.output_location);
    let reads_path = Path::new(&args.reads);

    // Create the output file so it errors if a bad output file is provided before computation
    let output_file = create_output_file(output_loc_path, "skim.r2f");

    // Create a mutex over a writer to allow multiple threads to write to the output file
    let output_writer = Mutex::new(BufWriter::new(output_file));

    let stats = Mutex::new((0, 0, 0.0, 0.0));

    info!("loading database at {:?}", database_path);
    let database = load_data_from_file::<Database>(database_path);

    info!("computing lookup table...");
    let lookup_table = database.compute_loookup_table(args.n_fixed as u64);

    info!(
        "classifying reads with cutoff threshold {}...",
        10.0_f64.powi((args.exp_cutoff).neg())
    );
    let read_iter = get_fastq_iter_of_file(reads_path);
    let start_time = Instant::now();

    read_iter
        .par_bridge()
        .into_par_iter()
        .for_each(|record_result| match record_result {
            Err(_) => {
                warn!("error encountered while reading fastq file");
                warn!("skipping the read that caused the error")
            }
            Ok(record) => {
                let (classification, (hit_lookup_time, prob_calc_time)) =
                    database.classify(record.seq(), cutoff_threshold, args.n_fixed, &lookup_table);

                {
                    let mut stats = stats.lock().unwrap();

                    stats.0 += 1;
                    stats.1 += record.seq().len();
                    stats.2 += hit_lookup_time;
                    stats.3 += prob_calc_time;
                }

                // Write classification result to output file
                let mut writer = output_writer.lock().unwrap();
                match classification {
                    Some((file, taxid)) => {
                        writer
                            .write(format!("C\t{}\t{}\t{}\n", record.id(), taxid, file).as_bytes())
                            .expect("could not write to output file");
                    }
                    None => {
                        writer
                            .write(format!("U\t{}\t0\t-\n", record.id()).as_bytes())
                            .expect("could not write to output file");
                    }
                };
            }
        });

    // Log throughput statisitcs of classification
    let classify_time = start_time.elapsed().as_secs_f64();
    let stats = stats.into_inner().unwrap();
    info!("classification took: {} s", classify_time);
    info!(
        "{} total reads classified ({} reads/s)",
        stats.0,
        stats.0 as f64 / classify_time
    );
    info!(
        "{} total bp classified ({} Mbp/s)",
        stats.1,
        (stats.1 as f64 / classify_time) / 1_000_000.0
    );
    debug!(
        "total thread time spent looking up kmer hits: {} s",
        stats.2
    );
    debug!(
        "total thread time spent calculating probabilities: {} s",
        stats.3
    );

    output_writer
        .into_inner()
        .expect("could not reclaim file writer at the end of execution")
        .flush()
        .expect("could not write to output file");

    info!("done!");
}

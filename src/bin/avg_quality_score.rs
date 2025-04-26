use clap::Parser;
use skim::io::create_output_file;
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::get_fastq_iter_of_file;
use std::io::{BufWriter, Write};
use std::path::Path;
use tracing::{info, warn};

/// Writes the average quality score for each read to the output
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// Fastq reads
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let reads_path = Path::new(&args.reads);

    let mut output_writer = BufWriter::new(create_output_file(output_loc_path, "r2avg_qscore"));

    let fastq_reads_iter = get_fastq_iter_of_file(reads_path);

    info!("looping through reads at {:?}", reads_path);

    let avg_qscore_per_read = fastq_reads_iter
        .into_iter()
        .filter_map(|record_result| match record_result {
            Ok(read) => {
                let total_qscores = read.qual().len() as f64;
                let sum_qscores = read
                    .qual()
                    .iter()
                    .map(|qscore| *qscore as usize)
                    .sum::<usize>() as f64;
                Some((read.id().to_string(), sum_qscores / total_qscores))
            }
            Err(e) => {
                warn!("error encountered: {:?} - skipping...", e);
                None
            }
        })
        .collect::<Vec<(String, f64)>>();

    let sum_avg_qscore = avg_qscore_per_read
        .iter()
        .map(|(_, qscore)| *qscore)
        .sum::<f64>();
    info!(
        "average quality score: {}",
        sum_avg_qscore / avg_qscore_per_read.len() as f64
    );

    for (readid, avg_qscore) in avg_qscore_per_read {
        output_writer
            .write(format!("{}\t{}\n", readid, avg_qscore).as_bytes())
            .expect("could not write to output file");
    }

    output_writer.flush().unwrap();

    info!("done!");
}

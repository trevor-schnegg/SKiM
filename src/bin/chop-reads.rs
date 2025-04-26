use bio::io::{fasta, fastq};
use clap::Parser;
use skim::io::create_output_file;
use skim::tracing::start_skim_tracing_subscriber;
use skim::utility::{get_fasta_iter_of_file, get_fastq_iter_of_file};
use std::path::Path;
use tracing::info;

/// Chop the input reads to the desired length
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, action)]
    /// Specify this flag if the file is a fasta file
    fasta: bool,

    #[arg(short, long, default_value_t = 180)]
    /// Maximum length of the read
    length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// Fasta reads file to chop
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_skim_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let chop_length = args.length;
    let reads_path = Path::new(&args.reads);

    if args.fasta {
        let output_file = create_output_file(output_loc_path, "chopped.fasta");
        let mut writer = fasta::Writer::new(output_file);

        let mut fasta_reads_iter = get_fasta_iter_of_file(reads_path);

        while let Some(Ok(read)) = fasta_reads_iter.next() {
            let seq = if read.seq().len() < chop_length {
                read.seq()
            } else {
                &read.seq()[..chop_length]
            };
            writer.write(read.id(), read.desc(), seq).unwrap();
        }
    } else {
        let output_file = create_output_file(output_loc_path, "chopped.fastq");
        let mut writer = fastq::Writer::new(output_file);

        let mut fastq_reads_iter = get_fastq_iter_of_file(reads_path);

        while let Some(Ok(read)) = fastq_reads_iter.next() {
            let (seq, qual) = if read.seq().len() < chop_length {
                (read.seq(), read.qual())
            } else {
                (&read.seq()[..chop_length], &read.qual()[..chop_length])
            };
            writer.write(read.id(), read.desc(), seq, qual).unwrap();
        }
    }

    info!("done!");
}

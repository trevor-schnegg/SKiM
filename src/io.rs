use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::any::type_name;
use std::fs::File;
use std::io::BufWriter;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tracing::{info, warn};

pub fn create_output_file(path: &Path, extension: &str) -> File {
    let file_path = if path.is_dir() {
        path.join(extension)
    } else {
        path.with_extension(extension)
    };

    info!("creating output file {:?}", file_path);

    File::create(file_path).expect("could not create output file")
}

pub fn split_string_to_taxid(line: String) -> Result<(String, usize), String> {
    let mut column_iter = line.split("\t");
    let file = column_iter.next().unwrap().to_string();

    // If there is no tab character, .next() will not exist
    match column_iter.next() {
        None => Err("line did not have a tab character".to_string()),
        Some(str) => {
            // Try to parse as a usize
            match str.parse::<usize>() {
                Ok(taxid) => Ok((file, taxid)),
                Err(e) => Err(e.to_string()),
            }
        }
    }
}

pub fn load_string2taxid(string2taxid: &Path) -> Vec<(String, usize)> {
    let file = File::open(string2taxid).expect(&*format!(
        "could not read string2taxid tsv at {:?}",
        string2taxid
    ));
    let reader = BufReader::new(file).lines();

    reader
        .enumerate()
        .filter_map(|(line_num, line)| {
            // Try to get the line from the file
            match line {
                Ok(line) => {
                    // Try to parse the line into a (file, taxid) tuple
                    match split_string_to_taxid(line) {
                        Ok(tuple) => Some(tuple),
                        Err(msg) => {
                            warn!(
                                "line {} from {:?} had an error :{}. skipping...",
                                line_num, string2taxid, msg
                            );
                            None
                        }
                    }
                }
                Err(e) => {
                    warn!(
                        "line {} from {:?} had an error :{}. skipping...",
                        line_num, string2taxid, e
                    );
                    None
                }
            }
        })
        .collect_vec()
}

// Takes a file (already opened) as an input
// All binaries open files at the start of execution, if needed.
// All such binaries should error early in execution if an improper path is provided.
pub fn dump_data_to_file<T: Serialize>(data: &T, file: File) -> bincode::Result<()> {
    let buf_writer = BufWriter::new(file);
    bincode::serialize_into(buf_writer, data)
}

// Takes a path (not opened) as an input
// All binaries that need to load data will do so at the start of execution.
// All such binaries will error here if an improper path is provided.
pub fn load_data_from_file<T: for<'a> Deserialize<'a>>(path: &Path) -> T {
    let buf_reader =
        BufReader::new(File::open(path).expect(&*format!("could not open file at {:?}", path)));
    bincode::deserialize_from(buf_reader).expect(&*format!(
        "failed to deserialize data at {:?} into {}",
        path,
        type_name::<T>()
    ))
}

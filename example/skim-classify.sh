#!/bin/bash

# Make sure that SKiM is compiled before running!

# General notes:
#	- For any binary (BIN) in this project, you can run `BIN -h` or `BIN --help` to see all options, arguments, and their descriptions.
#	- You can set the environment variable `RUST_LOG=debug` to print DEBUG messages to stdout

mkdir -p ./classification

# Perform classification on the test reads. This will create a readid2file (.r2f) file.
RUST_LOG=debug ../target/release/skim-classify -o ./classification/ ./database/skim.db ./test_reads/reads.fastq


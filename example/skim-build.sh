#!/bin/bash

# Make sure that SKiM is compiled before running!

# General notes:
#	- For any binary (BIN) in this project, you can run `BIN -h` or `BIN --help` to see all options, arguments, and their descriptions.
#	- You can set the environment variable `RUST_LOG=debug` to print DEBUG messages to stdout

mkdir -p ./database

# First create a file2taxid (.f2t) file from the reference files
../target/release/skim-file2taxid -a ./reference/accession2taxid -o ./database/ ./reference/sequences/

# Next, find the pairwise distances (.pd) for all reference files
RUST_LOG=debug ../target/release/skim-pairwise-distances -o ./database/ ./database/skim.f2t ./reference/sequences/

# Then, create an ordered file2taxid (.o.f2t) using the pairwise distances (.pd)
RUST_LOG=debug ../target/release/skim-order -o ./database/ ./database/skim.pd

# Finally, create a skim database (.db) from the ordered file2taxid (.o.f2t)
RUST_LOG=debug ../target/release/skim-build -o ./database/ ./database/skim.o.f2t ./reference/sequences/

# Optionally, create a lossy compressed database (.cdb) from the database (.db)
# This is left commented out because it is not worth doing on such a small example.
# In fact, it may not do anything based on how lossy compression is performed.

# RUST_LOG=debug ../target/release/skim-build-lossy -o ./database/ 1 ./database/skim.db


# SKiM
SKiM (Short K-mers in Metagenomics) is a metagenomic classifier implemented in Rust. The key feature of SKiM is that it uses short $k$-mers ($k=15$ or $k=16$) to perform classification. Short $k$-mers, in combination with data compression techniques and statistical correction, allow SKiM to perform fast and accurate metagenomic classification over large collections of assemblies with a very small memory footprint. Originally, SKiM was designed for classifcation of ONT reads (paticularly during adaptive sampling), where efficient computational methods are critical and the relatively high error rate of sequencing is likely to impact long $k$-mers ($k\geq 31$) more than short ones. However, there is no reason to believe that SKiM would not work under other contexts.


## Requirements

 * Rust/Cargo version 1.79.0 or higher (see [rust-lang.org](www.rust-lang.org/tools/install))
 * Currently, SKiM has only been tested on Linux - although it may work for other operating systems with Rust appropriately installed


## Installation

First, clone the repository into your current directory and `cd` into the root of the repository:

```
https://github.com/trevor-schnegg/SKiM.git
cd SKiM/
```

Then, use `compile.sh` to compile the project using `cargo` (make sure an appropritate version of `cargo` is available):

```
./compile.sh
```

All SKiM binaries will be located at `target/release/skim-*` from the root of the repository.


## Usage

Each SKiM binary has help documentation that details all available options and their descriptions, found by running `skim-* -h`. Please refer to this help documentation to get a better understanding of what each option does.

Optionally, the environment variable `RUST_LOG=debug` can be set before running any binary to log otherwise hidden information about in-depth statistics of index construction and/or classification.

If using a pre-built database, skip to [classification](#classification).

### Index Construction

1. Create a file2taxid (.f2t) file for the input FASTA files. As an example:

    ```
    skim-file2taxid -o example -s seqid2taxid ref/
    ```

    This will look for FASTA files in the `ref/` directory, use `seqid2taxid` (a file of the form `<seqid>\t<taxid>`) to annotate FASTA files with tax ids, and write the file2taxid to `example.skim.f2t`. A few additional notes about this step:

    * A `seqid2taxid` file is **not** required (but see [here](#getting-a-seqid2taxid) for how to get one). If not provided, all tax ids will be set to 0, but SKiM will still be able to report the file that a read hits to (if any). In this way, SKiM can be used independent of a taxonomy.
    * Because SKiM uses short $k$-mers, FASTA files in `ref/` whose total length is too long may need to be split to ensure classification accuracy. If needed, SKiM will create the directory `ref/skim/` to store split FASTA files. If a FASTA file is too large, SKiM will try to split the file by sequence. If a sequence is still too large, it will then chop the sequence into fragments with overlap (adjustable with option `-l`).

2. Compute the pairwise distances (.pd) matrix from the file2taxid (.f2t) As an example:

    ```
    skim-pairwise-distances -o example example.skim.f2t ref/
    ```

    This takes the file2taxid `example.skim.f2t`, along with the original FASTA file directory `ref/`, and outputs the pairwise distance matrix to `example.skim.pd`. This can be a **very** computationally expensive step, requiring a significant amount of RAM and time.

3. Create an ordered file2taxid (.o.f2t) from a pairwise distance (.pd) matrix. The ordering is chosen to significantly reduce the resulting database size (and increase speed). As an example:

    ```
    skim-order -o example example.skim.pd
    ```

    This takes the pairwise distance file `example.skim.pd` outputs the ordered file2taxid to `example.skim.o.f2t`.

4. Finally, build a SKiM database (.db) from the ordered file2taxid (.o.f2t). As an example:

    ```
    skim-build -o example example.skim.o.f2t ref/
    ```

    This takes the ordered file2taxid `example.skim.o.f2t`, along with the original FASTA file directory `ref/`, and outputs the SKiM database to `example.skim.db`. This can also be a computationally expensive step, both in terms of RAM and time. The resulting database file (`example.skim.db` in this case) is the only file needed to perform classification and can be run anywhere that SKiM is installed.

#### Index Construction Example

A fully functional index construction example is provided in the `example/` directory from the root of the repository. The example can be run by entering this directory (`cd example/`) and running:

```
./skim-build.sh
```

The output for each binary is placed in `database/` within the `example/` directory.

#### Getting a Seqid2taxid

If using the NCBI taxonomy, the seqid2taxid files are known as accession2taxid files and can be downloaded from NCBI's [ftp site](ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/) (specifically the `nucl_*.accession2taxid.*` files). However, these files are likely to contain a significant number of accessions that are not in your custom database. **Still need to fill what to do about this.**

#### Other Construction Methods

In reference to the provided [index construction steps](#index-construction), the only inputs technically required for running step 4 is **any** file2taxid (.f2t) and its corresponding FASTA file directory. This introduces two options:

1. Skip steps 2 and 3 in index construction entirely. In this case, the resulting database will be **far** from its minimum possible size and from its maximum possible throughput. However, will perform the same classification-wise and the total index preparation time and resources would be reduced.

2. Skip steps 2 and 3, but create a custom ordering for the file2taxid. For example, assuming that a seqid2taxid was provided in step 1, the file2taxid could be sorted by taxid using the linux `sort` command. In theory, at least the FASTA files with the same tax id will end up next to each other in the file. This benefits the way that SKiM performs compression, and it is likely to provide at least some benefits over the semi-random ordering obtained from step 1. Again, the database is still likely to be far from its minimum possbible size and maximum possible throughput. But, it will perform the same classification-wise and the total index preparation time and resources would be reduced.

#### Modifying $k$-mer Size and Sub-sampling

If you plan on modifying $k$-mer size and/or sub-sampling parameters ($s$ and $t$), please make sure to provide the same options to all binaries in the index construction process. If the desired parameters are absent from one or more steps, the overall database may not be as optimized as it could be.

### Classification

Once a database (.db) file is obtained, classification produces a readid2file (.r2f) mapping. As an example:

```
skim-classify -e 9 -o example example.skim.db reads.fastq
```

This performs classification on the FASTQ reads `reads.fastq` using the database `example.skim.db`, a statistical cutoff threshold of $10^{-9}$, and writes the output to example.skim.r2f. A few additional notes:

* The `-e` parameter is used in the equation $10^{-e}$ and, as outlined in our paper, is the statistical significance required to consider a read classified. By default, $e=9$. However, if you'd like to increase the precision of classification (at the cost of some recall), you should increase this parameter to $e=12$, $e=15$, or even $e=18$.
* Although the output file is called a readid2file (.r2f), the output format follows Kraken2's output as closely as possible. Specifically, it is a tab-delimited file where the columns are (from right to left):
    1. `U` for unclassified or `C` for classified.
    2. The read id, from the FASTQ header.
    3. The assigned tax id (**warning**: this will be 0 if unclassified **or** a seqid2taxid was not provided when constructing the database).
    4. The file to which the read is classified to (or `-` if unclassified).

#### Classification Example

A fully functional classification example is provided in the `example/` directory from the root of the repository. The example can be run by entering this directory (`cd example/`), following [these instructions](#index-construction-example) to create the database, and then running:

```
./skim-classify.sh
```

The output of classification is placed in `classification/` within the `example/` directory.

## License

## How to Cite

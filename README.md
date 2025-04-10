# SKiM
SKiM (Short K-mers in Metagenomics) is a metagenomic classifier implemented in Rust. The key feature of SKiM is that it uses short $k$-mers ($k=15$ or $k=16$) to perform classification. Short $k$-mers, in combination with data compression techniques and statistical correction, allow SKiM to perform fast and accurate metagenomic classification over large collections of assemblies with a very small memory footprint. Originally, SKiM was designed for classifcation of ONT reads (paticularly during adaptive sampling), where efficient computational methods are critical and the relatively high error rate of sequencing is likely to impact long $k$-mers ($k\geq 31$) more than short ones. However, there is no reason to believe that SKiM would not work under other contexts.

## Requirements

## Installation

## Usage

## License

## How to Cite

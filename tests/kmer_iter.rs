use itertools::Itertools;
use skim::kmer_iter::KmerIter;

#[test]
fn non_canonical() {
    let sequence = "CGATTAAAGATAGAAATACACGNTGCGAGCAATCAAATT";
    // according to A = 00, C = 01, G = 10, T = 11 and a kmer size of 14
    let sequence_kmers: Vec<usize> = vec![
        0b_01_10_00_11_11_00_00_00_10_00_11_00_10_00,
        0b_10_00_11_11_00_00_00_10_00_11_00_10_00_00,
        0b_00_11_11_00_00_00_10_00_11_00_10_00_00_00,
        0b_11_11_00_00_00_10_00_11_00_10_00_00_00_11,
        0b_11_00_00_00_10_00_11_00_10_00_00_00_11_00,
        0b_00_00_00_10_00_11_00_10_00_00_00_11_00_01,
        0b_00_00_10_00_11_00_10_00_00_00_11_00_01_00,
        0b_00_10_00_11_00_10_00_00_00_11_00_01_00_01,
        0b_10_00_11_00_10_00_00_00_11_00_01_00_01_10,
        0b_11_10_01_10_00_10_01_00_00_11_01_00_00_00,
        0b_10_01_10_00_10_01_00_00_11_01_00_00_00_11,
        0b_01_10_00_10_01_00_00_11_01_00_00_00_11_11,
    ];

    assert_eq!(
        sequence_kmers,
        KmerIter::from(sequence.as_bytes(), 14, false, 14).collect_vec()
    );
}

#[test]
fn canonical() {
    let sequence = "CGATTAAAGATAGAAATACACGNTGCGAGCAATCAAATT";
    // according to A = 00, C = 01, G = 10, T = 11 and a kmer size of 14
    let sequence_kmers: Vec<usize> = vec![
        0b_01_10_00_11_11_00_00_00_10_00_11_00_10_00,
        0b_10_00_11_11_00_00_00_10_00_11_00_10_00_00,
        0b_00_11_11_00_00_00_10_00_11_00_10_00_00_00,
        0b_00_11_11_11_01_11_00_11_01_11_11_11_00_00,
        0b_11_00_00_00_10_00_11_00_10_00_00_00_11_00,
        0b_00_00_00_10_00_11_00_10_00_00_00_11_00_01,
        0b_00_00_10_00_11_00_10_00_00_00_11_00_01_00,
        0b_00_10_00_11_00_10_00_00_00_11_00_01_00_01,
        0b_01_10_11_10_11_00_11_11_11_01_11_00_11_01,
        0b_11_10_01_10_00_10_01_00_00_11_01_00_00_00,
        0b_00_11_11_11_10_00_11_11_10_01_11_01_10_01,
        0b_00_00_11_11_11_10_00_11_11_10_01_11_01_10,
    ];

    assert_eq!(
        sequence_kmers,
        KmerIter::from(sequence.as_bytes(), 14, true, 14).collect_vec()
    );
}

#[test]
fn minimizers() {
    let sequence = "CGATTAAAGATAGAAATACACGNTGCGAGCAATCAAATT";
    // according to A = 00, C = 01, G = 10, T = 11 and a kmer size of 14
    let sequence_kmers: Vec<usize> = vec![
        0b_00_11_11_00_00_00_10_00_11_00_10_00_00_00,
        0b_00_00_00_10_00_11_00_10_00_00_00_11_00_01,
        0b_00_00_10_00_11_00_10_00_00_00_11_00_01_00,
        0b_01_10_00_10_01_00_00_11_01_00_00_00_11_11,
    ];

    assert_eq!(
        sequence_kmers,
        KmerIter::from(sequence.as_bytes(), 16, false, 14).collect_vec()
    );
}

#[test]
fn canonical_plus_minimizers() {
    let sequence = "CGATTAAAGATAGAAATACACGNTGCGAGCAATCAAATT";
    // according to A = 00, C = 01, G = 10, T = 11 and a kmer size of 14
    let sequence_kmers: Vec<usize> = vec![
        0b_00_11_11_00_00_00_10_00_11_00_10_00_00_00,
        0b_00_11_11_11_01_11_00_11_01_11_11_11_00_00,
        0b_00_11_11_00_00_00_10_00_11_00_10_00_00_00,
        0b_00_11_11_11_01_11_00_11_01_11_11_11_00_00,
        0b_00_00_00_10_00_11_00_10_00_00_00_11_00_01,
        0b_00_00_10_00_11_00_10_00_00_00_11_00_01_00,
        0b_00_00_11_11_11_10_00_11_11_10_01_11_01_10,
    ];

    assert_eq!(
        sequence_kmers,
        KmerIter::from(sequence.as_bytes(), 16, true, 14).collect_vec()
    );
}

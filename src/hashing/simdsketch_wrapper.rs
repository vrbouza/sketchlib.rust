// Wrapper around the simd-sketch crate to provide sketching for DNA sequences

use crate::sketch::Sketch;
use needletail::{parse_fastx_file, parser::Format};
use packed_seq::PackedNSeqVec;
use simd_sketch::SketchParams;

/// Wrapper around simd_sketch to make DNA sketches
pub fn sketch_with_simd(
    name: &String,
    fastx1: &String,
    fastx2: &Option<String>,
    min_qual: u8, //est_coverage : usize,
    mut sketchers: Vec<SketchParams>,
) -> Sketch {
    let mut reader_peek =
        parse_fastx_file(fastx1).unwrap_or_else(|_| panic!("Invalid path/file: {}", fastx1));
    let seq_peek = reader_peek
        .next()
        .expect("Invalid FASTA/Q record")
        .expect("Invalid FASTA/Q record");
    let mut reads = false;
    if seq_peek.format() == Format::Fastq {
        reads = true;
        for is in sketchers.iter_mut() {
            // is.coverage = est_coverage;
            is.coverage = 10000;
        }
    } else {
        for is in sketchers.iter_mut() {
            is.coverage = 1;
        }
    }

    let mut seqs = PackedNSeqVec::default();
    let mut reader =
        parse_fastx_file(fastx1).unwrap_or_else(|_| panic!("Invalid path/file: {fastx1}"));

    if min_qual == 0 || !reads {
        while let Some(record) = reader.next() {
            seqs.push_ascii(&record.expect("Invalid FASTA/Q record").seq());
            seqs.push_ascii(b"N");
        }

        if fastx2.is_some() {
            let mut reader = parse_fastx_file(fastx2.as_ref().unwrap())
                .unwrap_or_else(|_| panic!("Invalid path/file"));
            while let Some(record) = reader.next() {
                seqs.push_ascii(&record.expect("Invalid FASTA/Q record").seq());
                seqs.push_ascii(b"N");
            }
        }
    } else {
        let min_qual_usize = min_qual as usize;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            seqs.push_from_ascii_and_quality(
                &seqrec.seq(),
                &seqrec.qual().unwrap(),
                min_qual_usize,
            );
        }

        if fastx2.is_some() {
            let mut reader = parse_fastx_file(fastx2.as_ref().unwrap())
                .unwrap_or_else(|_| panic!("Invalid second path/file"));
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid FASTA/Q record");
                seqs.push_from_ascii_and_quality(
                    &seqrec.seq(),
                    &seqrec.qual().unwrap(),
                    min_qual_usize,
                );
            }
        }
    }

    // Run the sketching // NOTE TODO: there is no filtering for the reads, apart from quality! I.e., the difference on counts is zero!
    let simdsketches = sketchers
        .iter()
        .map(|is| is.build().sketch(seqs.as_slice()))
        .collect::<Vec<_>>();

    // Get sketchlib.rust sketch objects
    Sketch::from_sketch_simd(&simdsketches, &name.to_string(), reads)
}

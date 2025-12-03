//! Wrapper around the simd-sketch crate to provide sketching for DNA sequences

use crate::sketch::Sketch;
use needletail::{parse_fastx_file, parser::Format};
use packed_seq::PackedNSeqVec;
use simd_sketch::SketchParams;
use std::path::Path;

/// Wrapper around simd_sketch to make DNA sketches
pub fn sketch_with_simd(
    name: &String,
    fastx1: &String,
    fastx2: &Option<String>,
    min_qual: u8,
    est_coverage: usize,
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
            is.coverage = est_coverage;
        }
    } else {
        for is in sketchers.iter_mut() {
            is.coverage = 1;
        }
    }

    let seqs: PackedNSeqVec;
    let seqs2: Option<PackedNSeqVec>;

    if min_qual == 0 || !reads {
        seqs = PackedNSeqVec::from_fastx(Path::new(fastx1));

        if fastx2.is_some() {
            seqs2 = Some(PackedNSeqVec::from_fastx(Path::new(fastx2.as_ref().unwrap())));
        } else {
            seqs2 = None;
        }
    } else {
        seqs = PackedNSeqVec::from_fastq_with_quality(Path::new(fastx1), min_qual);

        if fastx2.is_some() {
            seqs2 = Some(PackedNSeqVec::from_fastq_with_quality(Path::new(fastx2.as_ref().unwrap()), min_qual));
        } else {
            seqs2 = None;
        }
    }
    
    // // Run the sketching
    let simdsketches = sketchers
        .iter()
        .map(|is| {
            if seqs2.is_some() {
                is.build().sketch_seqs(&[seqs.as_slice(), seqs2.as_ref().unwrap().as_slice()])
            } else {
                is.build().sketch(seqs.as_slice())
            }
        }).collect::<Vec<_>>();

    // Get sketchlib.rust sketch objects
    Sketch::from_sketch_simd(&simdsketches, &name.to_string(), reads)
}

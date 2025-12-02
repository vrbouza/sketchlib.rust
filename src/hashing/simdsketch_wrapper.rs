use simd_sketch::{Sketch as Sketch_simd, SketchParams};
use needletail::{parser::Format, parse_fastx_file};
use packed_seq::{PackedSeqVec, PackedNSeqVec, SeqVec};
use crate::sketch::Sketch;

/// Wrapper around simd_sketch to make DNA sketches
pub fn sketch_with_simd(name: &String, fastx1: &String, fastx2: &Option<String>, min_qual: u8, 
                        min_count: u16, mut sketchers: Vec<SketchParams>) -> Sketch {

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
            if min_count != 0 {
                is.duplicate = true;
            }
            is.coverage = 30;
        }
    } else {
        for is in sketchers.iter_mut() {
            is.duplicate = false;
        }
    }
    let simdsketches: Vec<Sketch_simd>;
    let mut reader =
        parse_fastx_file(fastx1).unwrap_or_else(|_| panic!("Invalid path/file: {fastx1}"));
    if min_qual == 0 || !reads {
        let mut seqs = vec![];
        while let Some(record) = reader.next() {
            seqs.push(PackedNSeqVec::from_ascii(&record.expect("Invalid FASTA/Q record").seq()));
        }

        if fastx2.is_some() {
            let mut reader =
                parse_fastx_file(fastx2.as_ref().unwrap()).unwrap_or_else(|_| panic!("Invalid path/file"));
            while let Some(record) = reader.next() {
                seqs.push(PackedNSeqVec::from_ascii(&record.expect("Invalid FASTA/Q record").seq()));
            }
        }

        let seqslices = seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>();

        // Run the sketching // NOTE TODO: there is no filtering for the reads, apart from quality! I.e., the difference on counts is zero!
        simdsketches = sketchers.iter().map(|is| {
            is.build().sketch_seqs(&seqslices)
        }).collect::<Vec<_>>();
    } else {
        let mut seqs = vec![];
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            let mut tmpseq = Vec::new();
            for (base, qual) in seqrec.seq().iter().zip(seqrec.qual().unwrap()) {
                if *qual >= min_qual {
                    tmpseq.push(*base);
                } else {
                    seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                    tmpseq.clear();
                }
            }
            if !tmpseq.is_empty() {
                seqs.push(PackedSeqVec::from_ascii(&tmpseq));
            }
        }

        if fastx2.is_some() {
            let mut reader =
                parse_fastx_file(fastx2.as_ref().unwrap()).unwrap_or_else(|_| panic!("Invalid second path/file"));
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid FASTA/Q record");
                let mut tmpseq = Vec::new();
                for (base, qual) in seqrec.seq().iter().zip(seqrec.qual().unwrap()) {
                    if *qual >= min_qual {
                        tmpseq.push(*base);
                    } else {
                        seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                        tmpseq.clear();
                    }
                }
                if !tmpseq.is_empty() {
                    seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                }
            }
        }
        let seqslices = seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>();

        // Run the sketching // NOTE TODO: there is no filtering for the reads, apart from quality! I.e., the difference on counts is zero!
        simdsketches = sketchers.iter().map(|is| {
            is.build().sketch_seqs(&seqslices)
        }).collect::<Vec<_>>();
    }

    // Get sketchlib.rust sketch objects
    Sketch::from_sketch_simd(&simdsketches, &name.to_string(), reads)
}
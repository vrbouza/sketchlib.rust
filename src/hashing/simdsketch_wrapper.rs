use simd_sketch::{Sketch as Sketch_simd, SketchParams};
use needletail::{parser::Format, parse_fastx_file};
use packed_seq::{PackedSeqVec, SeqVec};
use crate::hashing::{valid_base, encode_base};
use crate::sketch::Sketch;


pub fn sketch_with_simd(name: &String, fastx1: &String, fastx2: &Option<String>, min_qual: u8, 
                        min_count: u16, mut sketchers: Vec<SketchParams>) -> Sketch {

    // get k, ss, rc, from sketchers

    let mut reader_peek =
        parse_fastx_file(fastx1).unwrap_or_else(|_| panic!("Invalid path/file: {}", fastx1));
    let seq_peek = reader_peek
        .next()
        .expect("Invalid FASTA/Q record")
        .expect("Invalid FASTA/Q record");
    let mut reads = false;
    if seq_peek.format() == Format::Fastq && min_count != 0 {
        reads = true;
        for is in sketchers.iter_mut() {
            is.duplicate = true;
        }
    } else {
        for is in sketchers.iter_mut() {
            is.duplicate = false;
        }
    }

    let mut seqs = vec![];
    let mut reader =
        parse_fastx_file(fastx1).unwrap_or_else(|_| panic!("Invalid path/file: {fastx1}"));
    if min_qual == 0 {
        while let Some(record) = reader.next() {
            let mut tmpseq = Vec::new();
            let seqrec = record.expect("Invalid FASTA/Q record");
            for base in seqrec.seq().iter() {
                if valid_base(*base) {
                    let encoded_base = encode_base(*base);
                    tmpseq.push(encoded_base)
                } else {
                    seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                    tmpseq.clear();
                }
            }
            if !tmpseq.is_empty() {
                seqs.push(PackedSeqVec::from_ascii(&tmpseq));
            }
        } 
    } else {
        while let Some(record) = reader.next() {
            let mut tmpseq = Vec::new();
            let seqrec = record.expect("Invalid FASTA/Q record");
            if let Some(quals) = seqrec.qual() {
                for (base, qual) in seqrec.seq().iter().zip(quals) {
                    if *qual >= min_qual {
                        if valid_base(*base) {
                            let encoded_base = encode_base(*base);
                            tmpseq.push(encoded_base)
                        } else {
                            seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                            tmpseq.clear();
                        }
                    } else {
                        seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                        tmpseq.clear();
                    }
                }
            } else {
                for base in seqrec.seq().iter() {
                    if valid_base(*base) {
                        let encoded_base = encode_base(*base);
                        tmpseq.push(encoded_base)
                    } else {
                        seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                        tmpseq.clear();
                    }
                }
            }

            if !tmpseq.is_empty() {
                seqs.push(PackedSeqVec::from_ascii(&tmpseq));
            }
        } 
    }


    if fastx2.is_some() {
        let mut reader =
            parse_fastx_file(fastx2.as_ref().unwrap()).unwrap_or_else(|_| panic!("Invalid path/file: {fastx1}"));
        while let Some(record) = reader.next() {
            let mut tmpseq = Vec::new();
            let seqrec = record.expect("Invalid FASTA/Q record");
            if let Some(quals) = seqrec.qual() {
                for (base, qual) in seqrec.seq().iter().zip(quals) {
                    if *qual >= min_qual {
                        if valid_base(*base) {
                            let encoded_base = encode_base(*base);
                            tmpseq.push(encoded_base)
                        } else {
                            seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                            tmpseq.clear();
                        }
                    } else {
                        seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                        tmpseq.clear();
                    }
                }
            } else {
                for base in seqrec.seq().iter() {
                    if valid_base(*base) {
                        let encoded_base = encode_base(*base);
                        tmpseq.push(encoded_base)
                    } else {
                        seqs.push(PackedSeqVec::from_ascii(&tmpseq));
                        tmpseq.clear();
                    }
                }
            }

            if !tmpseq.is_empty() {
                seqs.push(PackedSeqVec::from_ascii(&tmpseq));
            }
        }
    }

    let seqslices = seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>();

    // Run the sketching // NOTE TODO: there is no filtering for the reads, apart from quality! I.e., the difference on counts is zero!
    let simdsketches: Vec<Sketch_simd> = sketchers.iter().map(|is| {
        is.build().sketch_seqs(&seqslices)
    }).collect::<Vec<_>>();

    // Get sketchlib.rust sketch objects
    Sketch::from_sketch_simd(&simdsketches, &name.to_string(), reads)
}
//! The class to support .ski creation, reading and writing, containing an inverted
//! index of multiple sketches
use std::fmt;
use std::sync::mpsc;

extern crate needletail;
use hashbrown::HashMap;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use rayon::prelude::*;
use roaring::{RoaringBitmap, RoaringTreemap};
use serde::{Deserialize, Serialize};

use super::hashing::{
    HashType, simdsketch_wrapper::sketch_with_simd
};
use crate::distances::distance_matrix::square_to_condensed;
use crate::io::InputFastx;
use crate::sketch::{sketch_datafile::SketchArrayWriter, BBITS};
use crate::utils::get_progress_bar;
use anyhow::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};

type InvSketches = (Vec<Vec<u16>>, Vec<String>);

use simd_sketch::{SketchParams, Sketch as Sketch_simd, SketchAlg, BitSketch};

/// An inverted index and associated metadata
#[derive(Serialize, Deserialize, Default, Clone, PartialEq)]
pub struct Inverted {
    index: Vec<HashMap<u16, RoaringBitmap>>,
    n_samples: usize,
    sample_names: Vec<String>,
    kmer_size: usize,
    sketch_version: String,
    rc: bool,
    hash_type: HashType,
}

impl Inverted {
    /// Sketch files without transposing bins, and invert the index. Files
    /// may be reordered via [`crate::io::reorder_input_files`]
    ///
    /// Optionally write untransposed samples to an skq file with `write_skq`.
    pub fn new(
        input_files: &[InputFastx],
        write_skq: Option<String>,
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        est_coverage: usize,
        quiet: bool,
    ) -> Self {
        log::info!("Creating sketches");
        let (sketches, names) = Self::sketch_files_inverted(
            input_files,
            file_order,
            k,
            sketch_size,
            seq_type,
            rc,
            min_count,
            min_qual,
            est_coverage,
            quiet,
        );
        if let Some(skq_file) = write_skq {
            log::info!("Writing bins for use with precluster as {skq_file}");
            let (bin_stride, kmer_stride, sample_stride) = (1, 1, sketch_size as usize);
            let mut skq_writer =
                SketchArrayWriter::new(&skq_file, bin_stride, kmer_stride, sample_stride);
            for sketch in &sketches {
                skq_writer.write_sketch(sketch);
            }
        }
        log::info!("Inverting sketch order");
        Self {
            index: Self::build_inverted_index(&sketches, sketch_size),
            n_samples: names.len(),
            sample_names: names,
            kmer_size: k,
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
            rc,
            hash_type: seq_type.clone(),
        }
    }

    /// Sketch query files in a compatible manner with the index,
    /// used when querying the index on the fly
    pub fn sketch_queries(
        &self,
        input_files: &[InputFastx],
        min_count: u16,
        min_qual: u8,
        est_coverage: usize,
        quiet: bool,
    ) -> InvSketches {
        let file_order: Vec<usize> = (0..input_files.len()).collect();
        Self::sketch_files_inverted(
            input_files,
            &file_order,
            self.kmer_size,
            self.index.len() as u64,
            &self.hash_type,
            self.rc,
            min_count,
            min_qual,
            est_coverage,
            quiet,
        )
    }

    /// Sample names in the index
    pub fn sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    /// Number of samples in the index
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    /// Sample name at the given index
    pub fn sample_at(&self, idx: usize) -> &str {
        &self.sample_names[idx]
    }

    /// K-mer length used for the index
    pub fn kmer(&self) -> usize {
        self.kmer_size
    }

    /// Sketch size used for the index
    pub fn sketch_size(&self) -> usize {
        self.index.len()
    }

    /// Saves to `file_prefix.ski`, using MessagePack as the serialisation format
    pub fn save(&self, file_prefix: &str) -> Result<(), Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Saving inverted index to {filename}");
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        rmp_serde::encode::write(&mut compress_writer, self)?;
        Ok(())
    }

    /// Loads from `file_prefix.ski`, using MessagePack as the serialisation format
    // NB MessagePack rather the CBOR uses here because of
    // https://github.com/enarx/ciborium/issues/96
    pub fn load(file_prefix: &str) -> Result<Self, Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Loading inverted index from {filename}");
        let ski_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let ski_obj: Self = rmp_serde::decode::from_read(decompress_reader)?;
        Ok(ski_obj)
    }

    /// Query a single sample against the index, returning the number of matches
    /// of bins against all samples in the index
    pub fn query_against_inverted_index(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut match_counts = vec![0; self.sample_names.len()];

        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                for sample_idx in matching_samples {
                    match_counts[sample_idx as usize] += 1;
                }
            }
        }
        match_counts
    }

    /// Return indexes of samples where all sketch bins are the same as a query
    pub fn all_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        matching_bits.insert_range(0..self.sample_names.len() as u32);
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits &= matching_samples;
            } else {
                matching_bits.clear();
                break;
            }
        }

        matching_bits.iter().collect()
    }

    /// Return indexes of samples where at least one sketch bins is the same as a query
    pub fn any_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits |= matching_samples;
            }
        }

        matching_bits.iter().collect()
    }

    /// Get the pairwise comparison list used by prefiler
    pub fn any_shared_bin_list(&self, quiet: bool) -> RoaringTreemap {
        let percent = false;
        let progress_bar = get_progress_bar(self.index.len(), percent, quiet);
        self.index
            .iter()
            .progress_with(progress_bar)
            .map(|bin| {
                bin.par_values()
                    .map(|hash_pres| {
                        let mut pair_map_hash = RoaringTreemap::new();
                        let samples_together: Vec<u32> = hash_pres.iter().collect();
                        for (i, sample1_idx) in samples_together.iter().enumerate() {
                            for sample2_idx in samples_together.iter().skip(i + 1) {
                                pair_map_hash.insert(square_to_condensed(
                                    *sample1_idx as usize,
                                    *sample2_idx as usize,
                                    self.n_samples,
                                ) as u64);
                            }
                        }
                        pair_map_hash
                    })
                    .reduce(RoaringTreemap::new, |pair_map_hash_a, pair_map_hash_b| {
                        pair_map_hash_a | pair_map_hash_b
                    }) // Reduction from pair_map_hash produces pair_map_bin
            })
            .reduce(|pair_map_all, pair_map_bin| pair_map_all | pair_map_bin)
            .unwrap()
        // Reduction from pair_map_bin produces pair_map_all
    }

    fn sketch_files_inverted(
        input_files: &[InputFastx],
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        est_coverage: usize,
        quiet: bool,
    ) -> InvSketches {
        let (tx, rx) = mpsc::channel();

        let percent = false;
        let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
        if *seq_type != HashType::DNA {
            panic!("Inverted index only available to DNA sequences.");
        }
        let sketchers = Some(vec![SketchParams {
            alg: SketchAlg::Bucket,
            rc: rc,
            k : k,
            s : sketch_size as usize,
            b : BBITS as usize,
            seed : 0,
            count : min_count as usize,
            coverage: 1,
            filter_empty: true,
            filter_out_n: true,
        }]);

        rayon::scope(|s| {
            s.spawn(move |_| {
                input_files
                    .par_iter()
                    .zip(file_order)
                    .progress_with(progress_bar)
                    .map(|((_, fastx1, fastx2), genome_idx)| {
                        // This could be nicer. It is essentially Sketch::from_sketch_simd, but w/o a couple things
                        // This vec has a single element always
                        let (_reads, vec) = sketch_with_simd(fastx1, fastx2, min_qual, est_coverage, sketchers.clone().unwrap());

                        let inputsketchs;
                        match &vec[0] {
                            Sketch_simd::BottomSketch(_sketch) => panic!("We only support BucketSketch at this moment."),
                            Sketch_simd::BucketSketch(sketch) => inputsketchs = &sketch.buckets,
                        }

                        let signs;
                        match inputsketchs {
                            BitSketch::B16(thevec) => signs = thevec,
                            _ => panic!("Only supporting 16 bits as bin size at this moment.")
                        }

                        (*genome_idx, signs.clone())
                    })
                    .for_each_with(tx, |tx, result| {
                        let _ = tx.send(result);
                    });
            });
        });

        let mut sketch_results = vec![Vec::new(); input_files.len()];
        while let Ok((genome_idx, sketch)) = rx.recv() {
            sketch_results[genome_idx] = sketch;
        }

        // Sample names in the correct order
        // (clones names, but reference would be annoying here)
        let mut sample_names: Vec<String> = vec!["".to_string(); input_files.len()];
        file_order
            .iter()
            .zip(input_files)
            .for_each(|(idx, (name, _, _))| sample_names[*idx] = name.to_string());

        (sketch_results, sample_names)
    }

    fn build_inverted_index(
        genome_sketches: &[Vec<u16>],
        sketch_size: u64,
    ) -> Vec<HashMap<u16, RoaringBitmap>> {
        // initialize inverted index structure
        let mut inverted_index: Vec<HashMap<u16, RoaringBitmap>> =
            vec![HashMap::new(); sketch_size as usize];

        // process each sketch
        // this could be parallelised over sketch bin, but probably not worth it
        // (NB doing par_iter on the below doesn't quite work as multiple mut borrows
        // of inverted_index needed)
        for (genome_idx, genome_signs) in genome_sketches.iter().enumerate() {
            for (i, hash) in genome_signs.iter().enumerate() {
                // add current genome to the inverted index at the current position
                inverted_index[i]
                    .entry(*hash)
                    .and_modify(|genome_list| {
                        genome_list.insert(genome_idx as u32);
                    })
                    .or_insert_with(|| {
                        let mut rb = RoaringBitmap::new();
                        rb.insert(genome_idx as u32);
                        rb
                    });
            }
        }
        inverted_index
    }
}

impl fmt::Debug for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sketch_version={}\nsequence_type={:?}\nsketch_size={}\nn_samples={}\nkmer={}\nrc={}\ninverted=true\n",
            self.sketch_version,
            self.hash_type,
            self.index.len(),
            self.sample_names.len(),
            self.kmer_size,
            self.rc,
        )?;
        let mut sizes = Vec::new();
        for bin in self.index.iter() {
            sizes.push(bin.len());
        }
        write!(
            f,
            "max_hashes_per_bin={}\nmin_hashes_per_bin={}\navg_hashes_per_bin={}",
            sizes.iter().max().unwrap(),
            sizes.iter().min().unwrap(),
            sizes.iter().sum::<usize>() as f64 / sizes.len() as f64
        )
    }
}

impl fmt::Display for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Name")?;
        for sketch in &self.sample_names {
            writeln!(f, "{sketch}")?;
        }
        Ok(())
    }
}

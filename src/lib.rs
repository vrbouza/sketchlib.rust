//! DOCS
//!

#![warn(missing_docs)]
use std::io::Write;
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::sketch::sketch_files;

pub mod multisketch;
use crate::multisketch::MultiSketch;

pub mod sketch_datafile;

pub mod distances;
use crate::distances::{calc_col_idx, calc_row_idx, DistType, DistanceMatrix};

pub mod io;
use crate::io::{get_input_list, parse_kmers, read_subset_names, set_ostream};

pub mod bloom_filter;
pub mod hashing;

/// Default k-mer size for sketching
pub const DEFAULT_KMER: usize = 17;
/// Chunk size in parallel distance calculations
pub const CHUNK_SIZE: usize = 1000;

#[doc(hidden)]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
        // simple_logger::init_with_level(log::Level::Trace).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    eprintln!("sketchlib.rust: fast biological distances at multiple k-mer lengths");
    let start = Instant::now();
    match &args.command {
        Commands::Sketch {
            seq_files,
            file_list,
            output,
            k_vals,
            k_seq,
            mut sketch_size,
            single_strand,
            min_count,
            min_qual,
            threads,
        } => {
            // An extra thread is needed for the writer. This doesn't 'overuse' CPU
            check_threads(*threads + 1);

            // Read input
            log::info!("Getting input files");
            let input_files = get_input_list(file_list, seq_files);
            let kmers = parse_kmers(k_vals, k_seq);
            // TODO this is very clunky, better replace fastx type
            let mut names: Vec<String> = input_files.iter().map(|x| x.0.to_string()).collect();
            // Build, merge
            let rc = !*single_strand;
            // Set expected sketchsize
            sketch_size = sketch_size.div_ceil(u64::BITS as u64);

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; threads:{}",
                kmers,
                sketch_size * u64::BITS as u64,
                threads
            );
            let mut sketches = sketch_files(
                &output,
                &input_files,
                &kmers,
                sketch_size,
                rc,
                *min_count,
                *min_qual,
            );
            let sketch_vec = MultiSketch::new(&mut names, &mut sketches, sketch_size, &kmers);
            sketch_vec
                .save_metadata(output)
                .expect("Error saving metadata");
        }
        Commands::Dist {
            ref_db,
            query_db,
            output,
            subset,
            kmer,
            threads,
        } => {
            check_threads(*threads);
            let mut output_file = set_ostream(output);

            let ref_db_name = if ref_db.ends_with(".skm") {
                &ref_db[0..ref_db.len() - 4]
            } else {
                ref_db.as_str()
            };
            log::info!("Loading sketch metadata from {}.skm", ref_db_name);
            let mut references = MultiSketch::load(ref_db_name)
                .expect(&format!("Could not read sketch metadata from {ref_db}.skm"));
            log::info!("Read sketches:\n{references:?}");

            log::info!("Loading sketch data from {}.skd", ref_db_name);
            if let Some(subset_file) = subset {
                let subset_names = read_subset_names(subset_file);
                references.read_sketch_data_block(ref_db_name, &subset_names);
            } else {
                references.read_sketch_data(ref_db_name);
            }

            // Set type of distances to use
            let k_idx;
            let dist_type = if let Some(k) = kmer {
                k_idx = references.get_k_idx(*k);
                DistType::Jaccard(*k)
            } else {
                k_idx = None;
                DistType::CoreAcc
            };
            log::info!("{dist_type}");

            let bar_style =
                ProgressStyle::with_template("{percent}% {bar:80.cyan/blue} eta:{eta}").unwrap();
            match query_db {
                None => {
                    // Self mode
                    log::info!("Calculating all ref vs ref distances");
                    let mut distances = DistanceMatrix::new(&references, None, dist_type);
                    let par_chunk = CHUNK_SIZE * distances.n_dist_cols();
                    let n = references.number_samples_loaded();
                    distances
                        .dists_mut()
                        .par_chunks_mut(par_chunk)
                        .progress_with_style(bar_style)
                        .enumerate()
                        .for_each(|(chunk_idx, dist_slice)| {
                            // Get first i, j index for the chunk
                            let start_dist_idx = chunk_idx * CHUNK_SIZE;
                            let mut i = calc_row_idx(start_dist_idx, n);
                            let mut j = calc_col_idx(start_dist_idx, i, n);
                            for dist_idx in 0..CHUNK_SIZE {
                                if let Some(k) = k_idx {
                                    let dist = references.jaccard_dist(i, j, k);
                                    dist_slice[dist_idx] = dist;
                                } else {
                                    let dist = references.core_acc_dist(i, j);
                                    dist_slice[dist_idx * 2] = dist.0;
                                    dist_slice[dist_idx * 2 + 1] = dist.1;
                                }

                                // Move to next index in upper triangle
                                j += 1;
                                if j >= n {
                                    i += 1;
                                    j = i + 1;
                                    // End of all dists reached (final chunk)
                                    if i >= (n - 1) {
                                        break;
                                    }
                                }
                            }
                        });

                    log::info!("Writing out in long matrix form");
                    write!(output_file, "{distances}").expect("Error writing output distances");
                }
                Some(db2) => {
                    // TODO Ref v query mode
                    log::info!("Calculating all ref vs query distances");
                    todo!()
                }
            }
        }
        Commands::Info {
            skm_file,
            sample_info,
        } => {
            let ref_db_name = if skm_file.ends_with(".skm") {
                &skm_file[0..skm_file.len() - 4]
            } else {
                skm_file.as_str()
            };
            let sketches = MultiSketch::load(ref_db_name).expect(&format!(
                "Could not read sketch metadata from {ref_db_name}.skm"
            ));
            println!("{sketches:?}");
            if *sample_info {
                log::info!("Printing sample info");
                println!("{sketches}");
            }
        }
    }
    let end = Instant::now();

    eprintln!(
        "🧬🖋️ sketchlib done in {}s",
        end.duration_since(start).as_secs()
    );
    log::info!("Complete");
}

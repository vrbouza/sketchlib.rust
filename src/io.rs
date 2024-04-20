use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use regex::Regex;

pub type InputFastx = (String, String, Option<String>);

pub fn read_input_fastas(seq_files: &[String]) -> Vec<InputFastx> {
    let mut input_files = Vec::new();
    // matches the file name (no extension) in a full path
    let re_path = Regex::new(r"^.+/(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    // matches the file name (no extension) with no path
    let re_name = Regex::new(r"^(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    for file in seq_files {
        let caps = re_path.captures(file).or(re_name.captures(file));
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, file.to_string(), None));
    }
    input_files
}

pub fn parse_kmers(k_list: &Option<Vec<usize>>, k_seq: &Option<Vec<usize>>) -> Vec<usize> {
    if k_list.is_some() && k_seq.is_some() {
        panic!("Only one of --k-vals or --k-seq should be specified");
    }

    let mut kmers = if let Some(k) = k_list {
        k.clone().to_vec()
    } else if let Some(k) = k_seq {
        (k[0]..=k[1]).step_by(k[2]).collect()
    } else {
        panic!("Must specify --k-vals or --k-seq");
    };

    kmers.sort_unstable();
    if !kmers.iter().all(|&k| k >= 5) {
        panic!("K-mers must be >=5");
    }

    kmers
}

/// Set a buffered stream to write to.
///
/// Either a file (if [`Some`]) or stdout otherwise (if [`None`]).
pub fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
    let out_writer = match oprefix {
        Some(prefix) => {
            let path = Path::new(prefix);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };
    BufWriter::new(out_writer)
}

/// Obtain a list of input files and names from command line input.
///
/// If `file_list` is provided, read each line as `name\tseq1\tseq2`, where
/// `seq2` is optional, and if present the reverse fastqs. Otherwise, treat
/// as fasta.
///
/// If `seq_files` are provided use [`read_input_fastas`].
pub fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<InputFastx> {
    // Read input
    match file_list {
        Some(files) => {
            let mut input_files: Vec<InputFastx> = Vec::new();
            let f = File::open(files).expect("Unable to open file_list");
            let f = BufReader::new(f);
            for line in f.lines() {
                let line = line.expect("Unable to read line in file_list");
                let fields: Vec<&str> = line.split_whitespace().collect();
                // Should be 2 entries for fasta, 3 for fastq
                let second_file = match fields.len() {
                    0..=1 => {
                        panic!("Unable to parse line in file_list")
                    }
                    2 => None,
                    3 => Some(fields[2].to_string()),
                    _ => {
                        panic!("Unable to parse line in file_list")
                    }
                };
                input_files.push((fields[0].to_string(), fields[1].to_string(), second_file));
            }
            input_files
        }
        None => read_input_fastas(seq_files.as_ref().unwrap()),
    }
}

pub fn read_subset_names(subset_file: &str) -> Vec<String> {
    let f = File::open(subset_file).expect(&format!("Unable to open {subset_file}"));
    let f = BufReader::new(f);
    let mut subset_names = Vec::new();
    for line in f.lines() {
        subset_names.push(line.unwrap());
    }
    log::info!("Read {} names to subset", subset_names.len());
    subset_names
}

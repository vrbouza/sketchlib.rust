use std::cmp::Ordering;
use std::fmt;
use std::sync::mpsc;

extern crate needletail;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::hashing::{nthash_iterator::NtHashIterator, HashType, RollHash};
use crate::bloom_filter::KmerFilter;
use crate::hashing::aahash_iterator::AaHashIterator;
use crate::io::InputFastx;
use crate::sketch_datafile::SketchArrayFile;
extern crate linfa;
extern crate linfa_clustering;
use linfa::DatasetBase;
use linfa::traits::{Fit, FitWith, Predict};
use linfa_clustering::{KMeansParams, KMeans, IncrKMeansError};

extern crate ndarray;
use ndarray::{Axis, array, s};
use ndarray_rand::rand::SeedableRng;
extern crate rand_xoshiro;
use rand_xoshiro::Xoshiro256Plus;

/*

pub fn fit_kmeans() {
    let observations = DatasetBase::from(data.clone());
    // Let's configure and run our K-means algorithm
    // We use the builder pattern to specify the hyperparameters
    // `n_clusters` is the only mandatory parameter.
    // If you don't specify the others (e.g. `n_runs`, `tolerance`, `max_n_iterations`)
    // default values will be used.
    let model = KMeans::params_with_rng(n_clusters, rng.clone())
        .tolerance(1e-2)
        .fit(&observations)
        .expect("KMeans fitted");

    // Once we found our set of centroids, we can also assign new points to the nearest cluster
    let new_observation = DatasetBase::from(array![[-9., 20.5]]);
    // Predict returns the **index** of the nearest cluster
    let dataset = model.predict(new_observation);
    // We can retrieve the actual centroid of the closest cluster using `.centroids()`
    let closest_centroid = &model.centroids().index_axis(Axis(0), dataset.targets()[0]);
    assert_abs_diff_eq!(closest_centroid.to_owned(), &array![-10., 20.], epsilon = 1e-1);
}*/

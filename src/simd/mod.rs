mod intrinsics;

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

use log::debug;
use packed_seq::{u32x8, ChunkIt, PackedNSeq, PaddedIt, Seq};
use seq_hash::KmerHasher;

/// Use the classic rotate-by-1 for backwards compatibility.
type FwdNtHasher = seq_hash::NtHasher<false, 1>;
type RcNtHasher = seq_hash::NtHasher<true, 1>;

#[derive(bincode::Encode, bincode::Decode, Debug)]
pub enum Sketch {
    BottomSketch(BottomSketch),
    BucketSketch(BucketSketch),
}

fn compute_mash_distance(j: f32, k: usize) -> f32 {
    assert!(j >= 0.0, "Jaccard similarity {j} should not be negative");
    // See eq. 4 of mash paper.
    let mash_dist = -(2. * j / (1. + j)).ln() / k as f32;
    assert!(
        mash_dist >= 0.0,
        "Bad mash distance {mash_dist} for jaccard similarity {j}"
    );
    // NOTE: Mash distance can be >1 when jaccard similarity is close to 0.
    // assert!(
    //     mash_dist <= 1.0,
    //     "Bad mash distance {mash_dist} for jaccard similarity {j}"
    // );
    // Distance 0 is computed as -log(1) and becomes -0.0.
    // This maximum fixes that.
    mash_dist.max(0.0)
}

impl Sketch {
    pub fn to_params(&self) -> SketchParams {
        match self {
            Sketch::BottomSketch(sketch) => SketchParams {
                alg: SketchAlg::Bottom,
                rc: sketch.rc,
                k: sketch.k,
                s: sketch.bottom.len(),
                b: 0,
                filter_empty: false,
                filter_out_n: false, // FIXME
            },
            Sketch::BucketSketch(sketch) => SketchParams {
                alg: SketchAlg::Bucket,
                rc: sketch.rc,
                k: sketch.k,
                s: sketch.buckets.len(),
                b: sketch.b,
                filter_empty: false,
                filter_out_n: false, // FIXME
            },
        }
    }
    pub fn jaccard_similarity(&self, other: &Self) -> f32 {
        match (self, other) {
            (Sketch::BottomSketch(a), Sketch::BottomSketch(b)) => a.jaccard_similarity(b),
            (Sketch::BucketSketch(a), Sketch::BucketSketch(b)) => a.jaccard_similarity(b),
            _ => panic!("Sketches are of different types!"),
        }
    }
    pub fn mash_distance(&self, other: &Self) -> f32 {
        let j = self.jaccard_similarity(other);
        let k = match self {
            Sketch::BottomSketch(sketch) => sketch.k,
            Sketch::BucketSketch(sketch) => sketch.k,
        };
        compute_mash_distance(j, k)
    }
}

/// Store only the bottom b bits of each input value.
#[derive(bincode::Encode, bincode::Decode, Debug)]
pub enum BitSketch {
    B32(Vec<u32>),
    B16(Vec<u16>),
    B8(Vec<u8>),
    B1(Vec<u64>),
}

impl BitSketch {
    fn new(b: usize, vals: Vec<u32>) -> Self {
        match b {
            32 => BitSketch::B32(vals),
            16 => BitSketch::B16(vals.into_iter().map(|x| x as u16).collect()),
            8 => BitSketch::B8(vals.into_iter().map(|x| x as u8).collect()),
            1 => BitSketch::B1({
                assert_eq!(vals.len() % 64, 0);
                vals.chunks_exact(64)
                    .map(|xs| {
                        xs.iter()
                            .enumerate()
                            .fold(0u64, |bits, (i, x)| bits | (((x & 1) as u64) << i))
                    })
                    .collect()
            }),
            _ => panic!("Unsupported bit width. Must be 1 or 8 or 16 or 32."),
        }
    }

    fn len(&self) -> usize {
        match self {
            BitSketch::B32(v) => v.len(),
            BitSketch::B16(v) => v.len(),
            BitSketch::B8(v) => v.len(),
            BitSketch::B1(v) => 64 * v.len(),
        }
    }
}

/// A sketch containing the `s` smallest k-mer hashes.
#[derive(bincode::Encode, bincode::Decode, Debug)]
pub struct BottomSketch {
    rc: bool,
    k: usize,
    bottom: Vec<u32>,
}

impl BottomSketch {
    /// Compute the similarity between two `BottomSketch`es.
    pub fn jaccard_similarity(&self, other: &Self) -> f32 {
        assert_eq!(self.rc, other.rc);
        assert_eq!(self.k, other.k);
        let a = &self.bottom;
        let b = &other.bottom;
        assert_eq!(a.len(), b.len());
        let mut intersection_size = 0;
        let mut union_size = 0;
        let mut i = 0;
        let mut j = 0;
        while union_size < a.len() {
            intersection_size += (a[i] == b[j]) as usize;
            let di = (a[i] <= b[j]) as usize;
            let dj = (a[i] >= b[j]) as usize;
            i += di;
            j += dj;
            union_size += 1;
        }

        return intersection_size as f32 / a.len() as f32;
    }

    pub fn mash_distance(&self, other: &Self) -> f32 {
        let j = self.jaccard_similarity(other);
        compute_mash_distance(j, self.k)
    }
}

/// A sketch containing the smallest k-mer hash for each remainder mod `s`.
#[derive(bincode::Encode, bincode::Decode, Debug)]
pub struct BucketSketch {
    rc: bool,
    k: usize,
    b: usize,
    pub buckets: BitSketch,
    /// Bit-vector indicating empty buckets, so the similarity score can be adjusted accordingly.
    empty: Vec<u64>,
}

impl BucketSketch {
    /// Compute the similarity between two `BucketSketch`es.
    pub fn jaccard_similarity(&self, other: &Self) -> f32 {
        assert_eq!(self.rc, other.rc);
        assert_eq!(self.k, other.k);
        assert_eq!(self.b, other.b);
        let both_empty = self.both_empty(other);
        // if both_empty > 0 {
        // debug!("Both empty: {}", both_empty);
        // }
        match (&self.buckets, &other.buckets) {
            (BitSketch::B32(a), BitSketch::B32(b)) => Self::inner_similarity(a, b, both_empty),
            (BitSketch::B16(a), BitSketch::B16(b)) => Self::inner_similarity(a, b, both_empty),
            (BitSketch::B8(a), BitSketch::B8(b)) => Self::inner_similarity(a, b, both_empty),
            (BitSketch::B1(a), BitSketch::B1(b)) => Self::b1_similarity(a, b, both_empty),
            _ => panic!("Bit width mismatch"),
        }
    }

    pub fn mash_distance(&self, other: &Self) -> f32 {
        let j = self.jaccard_similarity(other);
        compute_mash_distance(j, self.k)
    }

    fn inner_similarity<T: Eq>(a: &Vec<T>, b: &Vec<T>, both_empty: usize) -> f32 {
        assert_eq!(a.len(), b.len());
        let f = 1.0
            - std::iter::zip(a, b)
                .map(|(a, b)| (a != b) as u32)
                .sum::<u32>() as f32
                / (a.len() - both_empty) as f32;
        // Correction for accidental matches.
        let bb = (1usize << (size_of::<T>() * 8)) as f32;

        // Correction for accidental matches.
        // Take a max with 0 to avoid correcting into a negative jaccard similarity
        // for uncorrelated sketches.
        (bb * f - 1.0).max(0.0) / (bb - 1.0)
    }

    fn b1_similarity(a: &Vec<u64>, b: &Vec<u64>, both_empty: usize) -> f32 {
        assert_eq!(a.len(), b.len());
        let f = 1.0
            - std::iter::zip(a, b)
                .map(|(a, b)| (*a ^ *b).count_ones())
                .sum::<u32>() as f32
                / (64 * a.len() - both_empty) as f32;

        // Correction for accidental matches.
        // Take a max with 0 to avoid correcting into a negative jaccard similarity
        // for uncorrelated sketches.
        (2. * f - 1.).max(0.0)
    }

    fn both_empty(&self, other: &Self) -> usize {
        std::iter::zip(&self.empty, &other.empty)
            .map(|(a, b)| (a & b).count_ones())
            .sum::<u32>() as usize
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SketchAlg {
    Bottom,
    Bucket,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct SketchParams {
    /// Sketch algorithm to use. Defaults to bucket because of its much faster comparisons.
    pub alg: SketchAlg,
    /// When set, use forward instead of canonical k-mer hashes.
    pub rc: bool,
    /// k-mer size.
    pub k: usize,
    /// Bottom-s sketch, or number of buckets.
    pub s: usize,
    /// For bucket-sketch, store only the lower b bits.
    pub b: usize,
    /// For bucket-sketch, store a bitmask of empty buckets, to increase accuracy on small genomes.
    pub filter_empty: bool,

    pub filter_out_n: bool,
}

/// An object containing the sketch parameters.
///
/// Contains internal state to optimize the implementation when sketching multiple similar sequences.
pub struct Sketcher {
    params: SketchParams,
    rc_hasher: RcNtHasher,
    fwd_hasher: FwdNtHasher,
    factor: AtomicU64,
}

impl SketchParams {
    pub fn build(&self) -> Sketcher {
        let mut params = *self;
        // factor is pre-multiplied by 10 for a bit more fine-grained resolution.
        let factor;
        match params.alg {
            SketchAlg::Bottom => {
                // Clear out redundant value.
                params.b = 0;
                factor = 13;
            }
            SketchAlg::Bucket => {
                // To fill s buckets, we need ln(s)*s elements.
                // lg(s) is already a bit larger.
                factor = params.s.ilog2() as u64 * 5;
            }
        }
        if params.alg == SketchAlg::Bottom {}
        Sketcher {
            params,
            rc_hasher: RcNtHasher::new(params.k),
            fwd_hasher: FwdNtHasher::new(params.k),
            factor: AtomicU64::new(factor),
        }
    }

    /// Default sketcher that very fast at comparisons, but 20% slower at sketching.
    /// Use for >= 50000 seqs, and safe default when input sequences are > 500'000 characters.
    ///
    /// When sequences are < 100'000 characters, inaccuracies may occur due to empty buckets.
    pub fn default(k: usize) -> Self {
        SketchParams {
            alg: SketchAlg::Bucket,
            rc: true,
            k,
            s: 32768,
            b: 1,
            filter_empty: true,
            filter_out_n: false,
        }
    }

    /// Default sketcher that is fast at sketching, but somewhat slower at comparisons.
    /// Use for <= 5000 seqs, or when input sequences are < 100'000 characters.
    pub fn default_fast_sketching(k: usize) -> Self {
        SketchParams {
            alg: SketchAlg::Bucket,
            rc: true,
            k,
            s: 8192,
            b: 8,
            filter_empty: false,
            filter_out_n: false,
        }
    }
}

impl Sketcher {
    pub fn params(&self) -> &SketchParams {
        &self.params
    }

    /// Sketch a single sequence.
    pub fn sketch(&self, seq: impl Sketchable) -> Sketch {
        self.sketch_seqs(&[seq])
    }

    /// Sketch multiple sequence (fasta records) into a single sketch.
    pub fn sketch_seqs<'s>(&self, seqs: &[impl Sketchable]) -> Sketch {
        match self.params.alg {
            SketchAlg::Bottom => Sketch::BottomSketch(self.bottom_sketch(seqs)),
            SketchAlg::Bucket => Sketch::BucketSketch(self.bucket_sketch(seqs)),
        }
    }

    fn num_kmers<'s>(&self, seqs: &[impl Sketchable]) -> usize {
        seqs.iter()
            .map(|seq| seq.len() - self.params.k + 1)
            .sum::<usize>()
    }

    /// Return the `s` smallest `u32` k-mer hashes.
    /// Prefer [`Sketcher::sketch`] instead, which is much faster and just as
    /// accurate when input sequences are not too short.
    fn bottom_sketch<'s>(&self, seqs: &[impl Sketchable]) -> BottomSketch {
        // Iterate all kmers and compute 32bit nthashes.
        let n = self.num_kmers(seqs);
        let mut out = vec![];
        loop {
            let target = u32::MAX as usize * self.params.s / n;
            let factor = self.factor.load(Relaxed);
            let bound = (target as u128 * factor as u128 / 10 as u128).min(u32::MAX as u128) as u32;

            self.collect_up_to_bound(seqs, bound, &mut out);

            if bound == u32::MAX || out.len() >= self.params.s {
                out.sort_unstable();
                let old_len = out.len();
                out.dedup();
                let new_len = out.len();
                debug!("Deduplicated from {old_len} to {new_len}");
                if bound == u32::MAX || out.len() >= self.params.s {
                    out.resize(self.params.s, u32::MAX);

                    return BottomSketch {
                        rc: self.params.rc,
                        k: self.params.k,
                        bottom: out,
                    };
                }
            }

            let new_factor = factor + factor.div_ceil(4);
            let prev = self.factor.fetch_max(new_factor, Relaxed);
            debug!(
                "Found only {:>10} of {:>10} ({:>6.3}%)) Increasing factor from {factor} to {new_factor} (was already {prev})",
                out.len(),
                self.params.s,
                out.len() as f32 / self.params.s as f32,
            );
        }
    }

    /// s-buckets sketch. Splits the hashes into `s` buckets and returns the smallest hash per bucket.
    /// Buckets are determined via the remainder mod `s`.
    fn bucket_sketch<'s>(&self, seqs: &[impl Sketchable]) -> BucketSketch {
        // Iterate all kmers and compute 32bit nthashes.
        let n = self.num_kmers(seqs);
        let mut out = vec![];
        let mut buckets = vec![u32::MAX; self.params.s];
        loop {
            let target = u32::MAX as usize * self.params.s / n;
            let factor = self.factor.load(Relaxed);
            let bound = (target as u128 * factor as u128 / 10 as u128).min(u32::MAX as u128) as u32;

            debug!(
                "n {n:>10} s {} target {target:>10} factor {factor:>3} bound {bound:>10} ({:>6.3}%)",
                self.params.s,
                bound as f32 / u32::MAX as f32 * 100.0,
            );

            self.collect_up_to_bound(seqs, bound, &mut out);

            let mut empty = 0;
            if bound == u32::MAX || out.len() >= self.params.s {
                let m = FM32::new(self.params.s as u32);
                for &hash in &out {
                    let bucket = m.fastmod(hash);
                    buckets[bucket] = buckets[bucket].min(hash);
                }
                for &x in &buckets {
                    if x == u32::MAX {
                        empty += 1;
                    }
                }
                if bound == u32::MAX || empty == 0 {
                    if empty > 0 {
                        debug!("Found {empty} empty buckets.");
                    }
                    let empty = if empty > 0 && self.params.filter_empty {
                        debug!("Found {empty} empty buckets. Storing bitmask.");
                        buckets
                            .chunks(64)
                            .map(|xs| {
                                xs.iter().enumerate().fold(0u64, |bits, (i, x)| {
                                    bits | (((*x == u32::MAX) as u64) << i)
                                })
                            })
                            .collect()
                    } else {
                        vec![]
                    };

                    return BucketSketch {
                        rc: self.params.rc,
                        k: self.params.k,
                        b: self.params.b,
                        empty,
                        buckets: BitSketch::new(
                            self.params.b,
                            buckets.into_iter().map(|x| m.fastdiv(x) as u32).collect(),
                        ),
                    };
                }
            }

            let new_factor = factor + factor.div_ceil(4);
            let prev = self.factor.fetch_max(new_factor, Relaxed);
            debug!(
                "Found only {:>10} of {:>10} ({:>6.3}%, {empty:>5} empty) Increasing factor from {factor} to {new_factor} (was already {prev})",
                out.len(),
                self.params.s,
                out.len() as f32 / self.params.s as f32 * 100.,
            );
        }
    }

    fn collect_up_to_bound<'s>(&self, seqs: &[impl Sketchable], bound: u32, out: &mut Vec<u32>) {
        out.clear();
        if self.params.rc {
            for &seq in seqs {
                let hashes = seq.hash_kmers(&self.rc_hasher);
                collect_impl(bound, hashes, out);
            }
        } else {
            for &seq in seqs {
                let hashes = seq.hash_kmers(&self.fwd_hasher);
                collect_impl(bound, hashes, out);
            }
        }
        debug!(
            "Collect up to {bound:>10}: {:>9} ({:>6.3}%)",
            out.len(),
            out.len() as f32 / self.num_kmers(seqs) as f32 * 100.0
        );
    }
}

fn collect_impl(bound: u32, hashes: PaddedIt<impl ChunkIt<u32x8>>, out: &mut Vec<u32>) {
    let simd_bound = u32x8::splat(bound);
    let mut write_idx = out.len();
    let lane_len = hashes.it.len();
    let mut idx = u32x8::from(std::array::from_fn(|i| (i * lane_len) as u32));
    let max_idx = (8 * lane_len - hashes.padding) as u32;
    let max_idx = u32x8::splat(max_idx);
    hashes.it.for_each(|hashes| {
        let mask = hashes.cmp_lt(simd_bound);
        let in_bounds = idx.cmp_lt(max_idx);
        if write_idx + 8 > out.capacity() {
            out.reserve(out.capacity() + 8);
        }
        unsafe { intrinsics::append_from_mask(hashes, mask & in_bounds, out, &mut write_idx) };
        idx += u32x8::ONE;
    });

    unsafe { out.set_len(write_idx) };
}

pub trait Sketchable: Copy {
    fn len(&self) -> usize;
    fn hash_kmers<H: KmerHasher>(self, hasher: &H) -> PaddedIt<impl ChunkIt<u32x8>>;
}
impl Sketchable for &[u8] {
    fn len(&self) -> usize {
        Seq::len(self)
    }
    fn hash_kmers<H: KmerHasher>(self, hasher: &H) -> PaddedIt<impl ChunkIt<u32x8>> {
        hasher.hash_kmers_simd(self, 1)
    }
}
impl Sketchable for packed_seq::AsciiSeq<'_> {
    fn len(&self) -> usize {
        Seq::len(self)
    }
    fn hash_kmers<H: KmerHasher>(self, hasher: &H) -> PaddedIt<impl ChunkIt<u32x8>> {
        hasher.hash_kmers_simd(self, 1)
    }
}
impl Sketchable for packed_seq::PackedSeq<'_> {
    fn len(&self) -> usize {
        Seq::len(self)
    }
    fn hash_kmers<H: KmerHasher>(self, hasher: &H) -> PaddedIt<impl ChunkIt<u32x8>> {
        hasher.hash_kmers_simd(self, 1)
    }
}
impl<'s> Sketchable for PackedNSeq<'s> {
    fn len(&self) -> usize {
        Seq::len(&self.seq)
    }
    fn hash_kmers<'h, H: KmerHasher>(
        self,
        hasher: &'h H,
    ) -> PaddedIt<impl ChunkIt<u32x8> + use<'s, 'h, H>> {
        hasher.hash_valid_kmers_simd(self, 1)
    }
}

/// FastMod32, using the low 32 bits of the hash.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
struct FM32 {
    d: u64,
    m: u64,
}
impl FM32 {
    #[inline(always)]
    fn new(d: u32) -> Self {
        Self {
            d: d as u64,
            m: u64::MAX / d as u64 + 1,
        }
    }
    #[inline(always)]
    fn fastmod(self, h: u32) -> usize {
        let lowbits = self.m.wrapping_mul(h as u64);
        ((lowbits as u128 * self.d as u128) >> 64) as usize
    }
    #[inline(always)]
    fn fastdiv(self, h: u32) -> usize {
        ((self.m as u128 * h as u128) >> 64) as u32 as usize
    }
}

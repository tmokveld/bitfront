use criterion::{BenchmarkGroup, measurement::WallTime};
use serde::Deserialize;
use std::fs;
use std::hint::black_box;
use std::io;
use std::path::Path;

pub(super) struct Pair {
    pub query: Vec<u8>,
    pub target: Vec<u8>,
    pub distance: usize,
    pub cutoff: Option<usize>,
}

pub(super) fn measure<T>(
    group: &mut BenchmarkGroup<'_, WallTime>, directory: Option<&Path>, engine: &str,
    pairs: &[Pair], mut align: impl FnMut(&[u8], &[u8], Option<usize>) -> T,
) -> io::Result<()> {
    if let Some(directory) = directory {
        // Don't report last run's estimate if Criterion fails to save this one.
        // Criterion's saved baseline is in base/, not new/.
        match fs::remove_file(directory.join(engine).join("new/estimates.json")) {
            Err(err) if err.kind() == io::ErrorKind::NotFound => (),
            result => result?,
        }
    }
    group.bench_function(engine, |b| b.iter(|| {
        for pair in pairs {
            drop(black_box(align(black_box(pair.query.as_slice()), black_box(pair.target.as_slice()), pair.cutoff)));
        }
    }));
    Ok(())
}

// Criterion exposes estimates on disk, not through BenchmarkGroup. Only read
// the pinned version's median: WallTime is in ns and an iteration is one panel.
pub(super) fn median_per_pair(directory: &Path, engine: &str, pairs: usize) -> io::Result<f64> {
    #[derive(Deserialize)]
    struct Estimates { median: Estimate }
    #[derive(Deserialize)]
    struct Estimate { point_estimate: f64 }

    let path = directory.join(engine).join("new/estimates.json");
    let read = || -> io::Result<Estimates> { Ok(serde_json::from_slice(&fs::read(&path)?)?) };
    let estimates = read().map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    Ok(estimates.median.point_estimate / pairs as f64)
}

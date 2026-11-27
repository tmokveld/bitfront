use criterion::{BenchmarkGroup, measurement::WallTime};
use serde::Deserialize;
use std::fs;
use std::hint::black_box;
use std::io;
use std::path::Path;

pub(super) type Pair = (Vec<u8>, Vec<u8>);

pub(super) fn measure(
    group: &mut BenchmarkGroup<'_, WallTime>, directory: Option<&Path>, engine: &str,
    pairs: &[Pair], mut score: impl FnMut(&[u8], &[u8]) -> Option<usize>,
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
        for (query, target) in pairs {
            let _ = black_box(score(black_box(query.as_slice()), black_box(target.as_slice())));
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

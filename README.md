# Bitfront

<p align="center">
  <img src="assets/bitfront.png" width="128" alt="Bitfront logo">
</p>

Bitfront computes exact global, prefix, suffix, and infix Levenshtein alignments on byte slices. Ask for a score alone, a score with an edit path, or a score with just the matching prefix, suffix, or infix range. All support an optional inclusive cutoff.

## Add to your project

Bitfront is not published on crates.io. Use the repository's Git clone URL in your `Cargo.toml`:

```toml
[dependencies]
bitfront = { git = "https://github.com/tmokveld/bitfront.git" }
```

Bitfront requires Rust 1.85 or newer. Cargo fetches the [rust-edlib backend](https://crates.io/crates/rust-edlib) from crates.io for Myers scoring, range selection, and traceback.

## Comparing a read with a reference

It takes two edits to turn `ACGTAC` into `ACTTACC`: substitute `T` for `G`, then append `C`. Use `distance` to get that score, or `distance_at_most` to ask whether the read is within an edit limit.

```rust
fn main() -> Result<(), bitfront::Error> {
    let read = b"ACGTAC";
    let reference = b"ACTTACC";
    assert_eq!(bitfront::distance(read, reference)?, 2);

    // Two edits are allowed at the cutoff of 2, but not at the cutoff of 1.
    assert_eq!(bitfront::distance_at_most(read, reference, 2)?, Some(2));
    assert_eq!(bitfront::distance_at_most(read, reference, 1)?, None);
    Ok(())
}
```

## Matching part of a reference

Pass the read first and the reference second. These modes consume the whole read but only a contiguous part of the reference, which may be empty. Prefix and suffix scoring pin one end; `infix_distance` leaves both flanks free.

```rust
fn main() -> Result<(), bitfront::Error> {
    let reference = b"ACGTTAGC";
    // Left end: ACCT differs from the prefix ACGT by one substitution.
    assert_eq!(bitfront::prefix_distance(b"ACCT", reference)?, 1);
    assert_eq!(bitfront::prefix_distance_at_most(b"ACCT", reference, 0)?, None);
    // Right end: TAGC matches the suffix exactly.
    assert_eq!(bitfront::suffix_distance(b"TAGC", reference)?, 0);
    assert_eq!(bitfront::suffix_distance_at_most(b"TAGC", reference, 0)?, Some(0));
    // Anywhere: GTT is an internal match, so even with a cutoff of 1 the score is 0.
    assert_eq!(bitfront::infix_distance_at_most(b"GTT", reference, 1)?, Some(0));
    // Global scoring charges for the four leading bases that suffix scoring ignores.
    assert_eq!(bitfront::distance(b"TAGC", reference)?, 4);
    Ok(())
}
```

All scoring functions return a `Result`. With an inclusive cutoff `k`, `Ok(Some(d))` gives the exact distance and `Ok(None)` proves it exceeds `k`. `Error::InputTooLong` means the required dimensions or buffers cannot be accommodated, not that the cutoff was exceeded.

If you need to know *where* the query matched, `prefix_locate`, `suffix_locate`, and `infix_locate` return the score and one half-open byte range in the original target.

```rust
fn main() -> Result<(), bitfront::Error> {
    let target = b"ACGTTAGC";
    assert_eq!(bitfront::prefix_locate(b"ACG", target)?.target, 0..3);
    assert_eq!(bitfront::suffix_locate(b"TAGC", target)?.target, 4..8);
    let hit = bitfront::infix_locate_at_most(b"GTT", target, 1)?.unwrap();
    assert_eq!((hit.score, hit.target), (0, 2..5));
    Ok(())
}
```

## Getting an edit path

Use `align`, `prefix_align`, `suffix_align`, or `infix_align` to get one optimal alignment. Each has an `_at_most` variant: a hit includes the exact score and a path, and `Ok(None)` proves the score exceeds the cutoff. The `distance` functions do not construct paths.

```rust
fn main() -> Result<(), bitfront::Error> {
    let hit = bitfront::prefix_align(b"ACGT", b"ATCGTzz")?;
    assert_eq!(hit.score, 1);
    assert_eq!(hit.target, 0..5); // Only ATCGT is consumed.
    assert_eq!(hit.runs[1], bitfront::Run { op: bitfront::Op::Delete, len: 1 });
    assert_eq!(hit.cigar(bitfront::CigarFormat::Extended).to_string(), "1=1D3=");
    Ok(())
}
```

## How it works

Bitfront trims matching bytes from both ends for global distance, or just the pinned end for prefix and suffix distance. Its [WFA2-lib](https://github.com/Marco-Sola-Lab/WFA2-lib)-inspired wavefront searches in order of edit count and skips runs of equal bytes. When edits are dense, Bitfront can switch to [rust-edlib's](https://github.com/tmokveld/rust-edlib) bit-parallel Myers scorer.

## Tests

Run `cargo test --locked` for the library tests and API examples.

## Benchmark results

The [benchmark runner](tools/align-bench/README.md) compared Bitfront with the published [rust-edlib 0.5.1](https://crates.io/crates/rust-edlib/0.5.1) and [WFA2-lib via rust-wfa2 at `83d1b29`](https://github.com/tmokveld/rust-wfa2/tree/83d1b29cd06f8a7f8051b316d95ab0c5de7af073) on an Apple M5 Pro (arm64 macOS, rustc 1.99.1). WFA2 reused a `MemoryHigh` aligner with heuristics disabled. The comparison covers global, prefix, and suffix alignment, not infix.

Both tables use the same four fixed-seed DNA pairs per workload. Times are **median nanoseconds per pair** over 30 Criterion samples. Ratios use unrounded times: above 1× favors Bitfront.

### Score only

These calls return the distance without constructing an edit path.

| Search | Template bp | Query bp | Target bp | Edit rate | Distance | Bitfront ns/pair | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 64 | 64 | 64 | 0.00% | 0 | 4 | 4 | 29 | 1.06× | 7.13× |
| global | 1000 | 1001–1006 | 1000 | 5.00% | 48–49 | 1437 | 5905 | 2307 | 4.11× | 1.61× |
| global | 10000 | 9996–10030 | 10000 | 10.00% | 938–966 | 135603 | 213423 | 565583 | 1.57× | 4.17× |
| global | 256 | 253–261 | 1024 | 5.00% | 765–771 | 1061 | 5975 | 127793 | 5.63× | 120.46× |
| prefix | 256 | 253–261 | 1024 | 5.00% | 12–13 | 193 | 1603 | 398 | 8.30× | 2.06× |
| suffix | 256 | 253–261 | 1024 | 5.00% | 12–13 | 204 | 1601 | 14057 | 7.83× | 68.78× |
| global ≤ 32 (4/4 hit) | 256 | 256 | 256 | 12.50% | 32 | 621 | 1369 | 1151 | 2.20× | 1.85× |
| global ≤ 31 (0/4 hit) | 256 | 256 | 256 | 12.50% | 32 | 587 | 1365 | 1110 | 2.32× | 1.89× |
| prefix ≤ 32 (4/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 931 | 1535 | 1679 | 1.65× | 1.80× |
| prefix ≤ 31 (0/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 890 | 1506 | 1639 | 1.69× | 1.84× |
| suffix ≤ 32 (4/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 974 | 1556 | 35189 | 1.60× | 36.15× |
| suffix ≤ 31 (0/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 927 | 1530 | 34231 | 1.65× | 36.94× |

### Score + path

Each run includes the score and a path. Timings include constructing and dropping that result. The runner checks scores and replays each path against its target span before timing.

| Search | Template bp | Query bp | Target bp | Edit rate | Distance | Bitfront ns/pair | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 64 | 64 | 64 | 0.00% | 0 | 15 | 16 | 58 | 1.03× | 3.76× |
| global | 1000 | 1001–1006 | 1000 | 5.00% | 48–49 | 2246 | 11166 | 2770 | 4.97× | 1.23× |
| global | 10000 | 9996–10030 | 10000 | 10.00% | 938–966 | 334824 | 459252 | 588059 | 1.37× | 1.76× |
| global | 256 | 253–261 | 1024 | 5.00% | 765–771 | 4014 | 3992 | 132346 | 0.99× | 32.97× |
| prefix | 256 | 253–261 | 1024 | 5.00% | 12–13 | 325 | 4864 | 1921 | 14.97× | 5.92× |
| suffix | 256 | 253–261 | 1024 | 5.00% | 12–13 | 325 | 4983 | 15775 | 15.32× | 48.50× |
| global ≤ 32 (4/4 hit) | 256 | 256 | 256 | 12.50% | 32 | 786 | 2848 | 1324 | 3.63× | 1.69× |
| global ≤ 31 (0/4 hit) | 256 | 256 | 256 | 12.50% | 32 | 532 | 1404 | 1030 | 2.64× | 1.94× |
| prefix ≤ 32 (4/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 1165 | 4454 | 3213 | 3.82× | 2.76× |
| prefix ≤ 31 (0/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 898 | 1641 | 1573 | 1.83× | 1.75× |
| suffix ≤ 32 (4/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 1123 | 4675 | 36627 | 4.16× | 32.61× |
| suffix ≤ 31 (0/4 hit) | 256 | 256 | 1024 | 12.50% | 32 | 903 | 1656 | 34683 | 1.83× | 38.43× |

## License

[MIT](LICENSE).
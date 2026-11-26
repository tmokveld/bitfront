# Bitfront

<p align="center">
  <img src="assets/bitfront.png" width="128" alt="Bitfront logo">
</p>

Bitfront computes exact global, prefix, and suffix Levenshtein scores on byte slices, each with an optional cutoff.

## Add to your project

Bitfront is not published on crates.io. Use the repository's Git clone URL in your `Cargo.toml`:

```toml
[dependencies]
bitfront = { git = "https://github.com/tmokveld/bitfront.git" }
```

The library requires Rust 1.85 or newer and depends on [`rust-edlib` 0.3.1](https://crates.io/crates/rust-edlib/0.3.1). Both are pure Rust.

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

## Matching one end of a reference

Pass the read first and the reference second. The whole read is scored, but reference bases beyond the chosen prefix or suffix are free.

```rust
fn main() -> Result<(), bitfront::Error> {
    let reference = b"ACGTTAGC";
    // Left end: ACCT differs from the prefix ACGT by one substitution.
    assert_eq!(bitfront::prefix_distance(b"ACCT", reference)?, 1);
    assert_eq!(bitfront::prefix_distance_at_most(b"ACCT", reference, 0)?, None);
    // Right end: TAGC matches the suffix exactly.
    assert_eq!(bitfront::suffix_distance(b"TAGC", reference)?, 0);
    assert_eq!(bitfront::suffix_distance_at_most(b"TAGC", reference, 0)?, Some(0));
    // Global scoring charges for the four leading bases that suffix scoring ignores.
    assert_eq!(bitfront::distance(b"TAGC", reference)?, 4);
    Ok(())
}
```

An empty prefix or suffix is allowed. Neither function searches for a match in the middle of the target.

All scoring functions return a `Result`. With an inclusive cutoff `k`, `Ok(Some(d))` gives the exact distance and `Ok(None)` proves it exceeds `k`. `Error::InputTooLong` means the required dimensions or buffers cannot be accommodated, not that the cutoff was exceeded.

## How it works

Bitfront trims matching bytes from both ends for global distance, or just the pinned end for prefix and suffix distance. Its wavefront search is inspired by [WFA2-lib](https://github.com/Marco-Sola-Lab/WFA2-lib): it searches in order of edit count and skips runs of equal bytes. That makes it fast when only a few edits separate the inputs. When edits are dense, the wavefront has more diagonals to explore. Bitfront can switch to Myers' bit-vector algorithm, which [Edlib](https://github.com/Martinsos/edlib) uses to update many edit-distance cells per machine word. 

## Tests

Run `cargo test --locked` for the library tests and API examples.

## Benchmark results

The baselines are [`rust-edlib` 0.3.1](https://crates.io/crates/rust-edlib/0.3.1), a pure-Rust implementation of [Edlib](https://github.com/Martinsos/edlib), and [`rust-wfa2` at `83d1b29`](https://github.com/tmokveld/rust-wfa2/tree/83d1b29cd06f8a7f8051b316d95ab0c5de7af073), the Rust bindings used to run [WFA2-lib](https://github.com/tmokveld/WFA2-lib) in score-only mode with `MemoryHigh` and heuristics disabled.

Timings are **median nanoseconds per pair** over five samples on four fixed-seed DNA pairs per row. All three engines agreed on every score and cutoff result in the table before timing. The edit rate is the fraction of edited sites in the source, the distance column reports the actual optimal score or range.

The last two columns divide each baseline's median time by Bitfront's before rounding. Above 1× favors Bitfront, below 1× favors the baseline.

| Search | Source bp | Target bp | Nominal edit rate | Exact distance | Bitfront ns/pair | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 64 | 64 | 0.00% | 0 | 5 | 260 | 39 | 47.41× | 7.06× |
| global | 64 | 64 | 1.00% | 1 | 6 | 197 | 51 | 33.52× | 8.60× |
| global | 64 | 64 | 5.00% | 3 | 13 | 197 | 80 | 15.15× | 6.20× |
| global | 64 | 64 | 20.00% | 12–13 | 135 | 192 | 310 | 1.42× | 2.30× |
| global | 256 | 256 | 0.00% | 0 | 11 | 1529 | 38 | 136.01× | 3.41× |
| global | 256 | 256 | 1.00% | 3 | 36 | 1551 | 86 | 42.52× | 2.36× |
| global | 256 | 256 | 5.00% | 12–13 | 190 | 1562 | 315 | 8.23× | 1.66× |
| global | 256 | 256 | 20.00% | 41–46 | 757 | 1569 | 1780 | 2.07× | 2.35× |
| global | 1000 | 1000 | 0.00% | 0 | 37 | 6265 | 82 | 167.07× | 2.18× |
| global | 1000 | 1000 | 1.00% | 10 | 166 | 6424 | 285 | 38.67× | 1.72× |
| global | 1000 | 1000 | 5.00% | 48–49 | 1515 | 6243 | 2261 | 4.12× | 1.49× |
| global | 1000 | 1000 | 20.00% | 176–186 | 9248 | 16055 | 23237 | 1.74× | 2.51× |
| global | 10000 | 10000 | 0.00% | 0 | 299 | 64232 | 502 | 214.95× | 1.68× |
| global | 10000 | 10000 | 0.10% | 10 | 510 | 65143 | 813 | 127.77× | 1.60× |
| global | 10000 | 10000 | 1.00% | 100 | 5361 | 109831 | 8705 | 20.49× | 1.62× |
| global | 10000 | 10000 | 5.00% | 480–490 | 95417 | 165542 | 155964 | 1.73× | 1.63× |
| global | 10000 | 10000 | 10.00% | 938–966 | 142747 | 229443 | 554948 | 1.61× | 3.89× |
| global | 100000 | 100000 | 0.00% | 0 | 2808 | 626760 | 5603 | 223.22× | 2.00× |
| global | 100000 | 100000 | 0.10% | 100 | 9009 | 1045875 | 15186 | 116.10× | 1.69× |
| global | 100000 | 100000 | 1.00% | 995–998 | 393255 | 2285281 | 728417 | 5.81× | 1.85× |
| global | 1000000 | 1000000 | 0.00% | 0 | 28570 | 6351562 | 55081 | 222.32× | 1.93× |
| global | 1000000 | 1000000 | 0.01% | 100 | 37641 | 10261740 | 68866 | 272.62× | 1.83× |
| global ≤ 32 (4/4 hit) | 256 | 256 | 12.50% | 32 | 613 | 1505 | 1134 | 2.46× | 1.85× |
| global ≤ 31 (0/4 hit) | 256 | 256 | 12.50% | 32 | 591 | 1482 | 1094 | 2.51× | 1.85× |
| global | 256 | 1024 | 5.00% | 765–771 | 1969 | 6951 | 132685 | 3.53× | 67.37× |
| prefix | 256 | 1024 | 5.00% | 12–13 | 204 | 1567 | 400 | 7.70× | 1.96× |
| suffix | 256 | 1024 | 5.00% | 12–13 | 208 | 1562 | 14189 | 7.50× | 68.11× |
| prefix | 10000 | 40000 | 1.00% | 100 | 8903 | 100129 | 13840 | 11.25× | 1.55× |
| suffix | 10000 | 40000 | 1.00% | 100 | 9289 | 103661 | 4098812 | 11.16× | 441.24× |

## License

[MIT](LICENSE).
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

Bitfront requires Rust 1.85 or newer and depends on [`rust-edlib` 0.4.0](https://crates.io/crates/rust-edlib/0.4.0) from crates.io. Both libraries are pure Rust.

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

The baselines are [`rust-edlib` 0.4.0](https://crates.io/crates/rust-edlib/0.4.0), a pure-Rust implementation of [Edlib](https://github.com/Martinsos/edlib), and [`rust-wfa2` at `83d1b29`](https://github.com/tmokveld/rust-wfa2/tree/83d1b29cd06f8a7f8051b316d95ab0c5de7af073), the Rust bindings used to run [WFA2-lib](https://github.com/tmokveld/WFA2-lib) in score-only mode with `MemoryHigh` and heuristics disabled.

Timings are **median nanoseconds per pair** over 30 Criterion samples on four fixed-seed DNA pairs per row. All three engines agreed on every score and cutoff result in the table before timing. The edit rate is the fraction of edited sites in the source, the distance column reports the actual optimal score or range.

The last two columns divide each baseline's median time by Bitfront's before rounding. Above 1× favors Bitfront, below 1× favors the baseline.

| Search | Source bp | Target bp | Nominal edit rate | Exact distance | Bitfront ns/pair | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 64 | 64 | 0.00% | 0 | 4 | 193 | 28 | 52.92× | 7.79× |
| global | 64 | 64 | 1.00% | 1 | 5 | 194 | 48 | 35.73× | 8.76× |
| global | 64 | 64 | 5.00% | 3 | 12 | 196 | 80 | 16.79× | 6.81× |
| global | 64 | 64 | 20.00% | 12–13 | 132 | 195 | 319 | 1.48× | 2.42× |
| global | 256 | 256 | 0.00% | 0 | 7 | 1450 | 37 | 216.08× | 5.47× |
| global | 256 | 256 | 1.00% | 3 | 33 | 1456 | 85 | 43.79× | 2.56× |
| global | 256 | 256 | 5.00% | 12–13 | 172 | 1468 | 319 | 8.55× | 1.85× |
| global | 256 | 256 | 20.00% | 41–46 | 637 | 1468 | 1785 | 2.30× | 2.80× |
| global | 1000 | 1000 | 0.00% | 0 | 20 | 5970 | 80 | 304.19× | 4.05× |
| global | 1000 | 1000 | 1.00% | 10 | 145 | 6049 | 284 | 41.76× | 1.96× |
| global | 1000 | 1000 | 5.00% | 48–49 | 1467 | 5922 | 2251 | 4.04× | 1.53× |
| global | 1000 | 1000 | 20.00% | 176–186 | 7502 | 14744 | 23526 | 1.97× | 3.14× |
| global | 10000 | 10000 | 0.00% | 0 | 170 | 60620 | 496 | 357.47× | 2.93× |
| global | 10000 | 10000 | 0.10% | 10 | 304 | 60918 | 783 | 200.59× | 2.58× |
| global | 10000 | 10000 | 1.00% | 100 | 5281 | 100041 | 8654 | 18.94× | 1.64× |
| global | 10000 | 10000 | 5.00% | 480–490 | 91628 | 154776 | 157150 | 1.69× | 1.72× |
| global | 10000 | 10000 | 10.00% | 938–966 | 134198 | 219823 | 554362 | 1.64× | 4.13× |
| global | 100000 | 100000 | 0.00% | 0 | 1635 | 636695 | 5646 | 389.38× | 3.45× |
| global | 100000 | 100000 | 0.10% | 100 | 7161 | 1015455 | 15163 | 141.80× | 2.12× |
| global | 100000 | 100000 | 1.00% | 995–998 | 393015 | 2113052 | 731378 | 5.38× | 1.86× |
| global | 1000000 | 1000000 | 0.00% | 0 | 17090 | 6163412 | 54956 | 360.65× | 3.22× |
| global | 1000000 | 1000000 | 0.01% | 100 | 25677 | 9963177 | 66256 | 388.03× | 2.58× |
| global ≤ 32 (4/4 hit) | 256 | 256 | 12.50% | 32 | 590 | 1463 | 1135 | 2.48× | 1.92× |
| global ≤ 31 (0/4 hit) | 256 | 256 | 12.50% | 32 | 561 | 1452 | 1096 | 2.59× | 1.95× |
| prefix ≤ 32 (4/4 hit) | 256 | 1024 | 12.50% | 32 | 936 | 1482 | 1706 | 1.58× | 1.82× |
| prefix ≤ 31 (0/4 hit) | 256 | 1024 | 12.50% | 32 | 902 | 1458 | 1580 | 1.62× | 1.75× |
| suffix ≤ 32 (4/4 hit) | 256 | 1024 | 12.50% | 32 | 901 | 1515 | 34562 | 1.68× | 38.35× |
| suffix ≤ 31 (0/4 hit) | 256 | 1024 | 12.50% | 32 | 864 | 1500 | 33783 | 1.74× | 39.11× |
| global | 256 | 1024 | 5.00% | 765–771 | 753 | 6452 | 128054 | 8.57× | 170.08× |
| prefix | 256 | 1024 | 5.00% | 12–13 | 206 | 1568 | 407 | 7.63× | 1.98× |
| suffix | 256 | 1024 | 5.00% | 12–13 | 210 | 1589 | 14105 | 7.56× | 67.10× |
| prefix | 10000 | 40000 | 1.00% | 100 | 8264 | 102550 | 13641 | 12.41× | 1.65× |
| suffix | 10000 | 40000 | 1.00% | 100 | 8708 | 104181 | 4096672 | 11.96× | 470.43× |

## License

[MIT](LICENSE).
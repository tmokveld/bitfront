# Running the alignment benchmark

The benchmark compares Bitfront with direct score-only calls to `rust-edlib` 0.4.0 from crates.io and with WFA2-lib. Both Rust scorers receive the same byte slices. Bounded global, prefix, and suffix rows call `distance_at_most`, `prefix_distance_at_most`, and `suffix_distance_at_most`, respectively. WFA2 uses the same inclusive cutoff and frees only the unpinned target end for prefix and suffix searches.

The [CI benchmark](../../.github/workflows/benchmark.yml) builds WFA2-lib with a C/C++ toolchain, Make, CMake 3.16 or newer, and `libclang`. Bitfront itself does not need these tools.

From the repository root, run the locked release benchmark:

```sh
RUSTUP_TOOLCHAIN=1.99.1 sh tools/align-bench/run.sh
```

By default, Criterion takes 30 samples per engine and workload. One iteration scores all four fixed-seed pairs; the table divides Criterion’s median time per iteration by four. `--smoke` uses one pair, a smaller panel, and 10 samples per engine. Use `--samples N` (or `--sample-size N`) to override either default; Criterion requires at least 10. This changes the number of timing samples, not the number of input pairs.

The script writes the Markdown table to stdout and Criterion’s progress and confidence intervals to stderr. Criterion saves its estimates and samples under `tools/align-bench/target/criterion` by default (`CRITERION_HOME` overrides this). For native Criterion output and filtering, run `cargo bench --locked --manifest-path tools/align-bench/Cargo.toml --bench align -- --help`. The [published results](../../README.md#benchmark-results) use the default panel.

On macOS, if `bindgen` cannot find Homebrew LLVM, set `LIBCLANG_PATH="$(brew --prefix llvm)/lib"`.

To time one workload, pass `--case MODE SOURCE_BP RATE_BP [CUTOFF]`. `MODE` is `global`, `prefix`, or `suffix`; `SOURCE_BP` is the source length in bases, and `RATE_BP` is the number of edits per 10,000 source bytes. Omit `CUTOFF` for an unbounded search. For example, this scores a 10,000 bp source at 1% edits against a 40,000 bp target:

```sh
sh tools/align-bench/run.sh --case prefix 10000 100 --samples 50
```

Bounded cases replace selected query bases with `N`, which never occurs in the target. Each replacement therefore costs exactly one edit. The standard and smoke panels include hits at `CUTOFF=32` and misses at `31` for all three searches. The prefix and suffix cases use a 256 bp query and a 1,024 bp target, with the extra target sequence at the free end.

## Compare with a previous Bitfront commit

Pass `--compare REV` first to add the specified revision as a fourth engine. `REV` is a commit, tag, or branch in this local Git repository. For example, compare the current working tree (including uncommitted changes) with the initial release:

```sh
RUSTUP_TOOLCHAIN=1.99.1 sh tools/align-bench/run.sh --compare 23f3d9c639a6a9950eabef29853dcdb2288c2180
```

Append `--smoke` or `--case MODE SOURCE_BP RATE_BP [CUTOFF]` to use a smaller panel. For example, `sh tools/align-bench/run.sh --compare HEAD --samples 50 --case suffix 256 1250 31` times a bounded suffix miss. The `previous / Bitfront` column divides the commit's median time by the working tree's: above 1× means the working tree is faster. Both versions score the same byte slices in the same release binary. The runner checks each result against the unbounded rust-edlib score for the same search mode before timing. Criterion warms up and measures one engine at a time.

This option requires Python 3 and Git. It resolves `REV` to a full commit ID and loads that source through a local Git dependency, without switching branches or modifying either tracked lockfile. Cargo resolves both versions together in a temporary project seeded with the benchmark lockfile; the chosen commit must have compatible dependencies and expose `distance`, `prefix_distance`, `suffix_distance`, and their `_at_most` variants.

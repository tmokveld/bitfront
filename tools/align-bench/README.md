# Running the alignment benchmark

The benchmark compares Bitfront, direct calls to [the same rust-edlib release from crates.io](https://crates.io/crates/rust-edlib), and WFA2-lib on identical byte slices. By default it reports **score** and **path** rows separately for global, prefix and suffix alignment. Use `--task score`, `--task path` or `--task both` to select.

Score rows request no traceback. Path rows time production and drop of the owned result: Bitfront and rust-edlib both return one optimal target span and coalesced runs (`align_runs` for direct rust-edlib); WFA2 returns a packed SAM CIGAR. The heap comparison uses the same calls. WFA2 reuses one single-threaded `MemoryHigh` aligner per case and task, with heuristics disabled; only path rows use alignment scope. Before timing, each path is replayed against the original query and reported target range, and its edit cost is checked against the unbounded distance. Ties need not choose the same path or span.

The [CI benchmark](../../.github/workflows/benchmark.yml) builds WFA2-lib with a C/C++ toolchain, Make, CMake 3.16 or newer, and `libclang`. Bitfront itself does not need these tools.

From the repository root, run the locked release benchmark:

```sh
RUSTUP_TOOLCHAIN=1.99.1 sh tools/align-bench/run.sh
```

By default, Criterion takes 30 samples per engine and workload. One iteration aligns all four fixed-seed pairs; the table divides Criterion’s median time per iteration by four. `--smoke` uses one pair, a smaller panel, and 10 samples per engine. Use `--samples N` (or `--sample-size N`) to override either default; Criterion requires at least 10.

The script writes the timing table to stdout first, followed by a **Score memory** table, a **Path memory** table, or both, according to `--task`. Criterion’s progress and confidence intervals go to stderr. Estimates and samples are saved under `tools/align-bench/target/criterion` by default (`CRITERION_HOME` overrides this). Path group IDs end in `-path`; their direct rust-edlib series is `rust-edlib-runs`, distinct from older unit-path measurements. Score IDs are unchanged. For timing only, with native Criterion filters, run `cargo bench --locked --manifest-path tools/align-bench/Cargo.toml --bench align -- --help`. The [README results](../../README.md#benchmark-results) show score-only and score-plus-path timings separately for the same workloads.

Memory is measured in a separate process on macOS and glibc Linux; the timing process never loads allocation hooks. **Peak B** is the most live heap for one engine across the same seeded pairs, not a sum or a median. It counts allocator-usable bytes (including size-class rounding) for Rust and native C allocations. WFA2 starts with a fresh aligner per row; its construction, retained buffers and owned results all count. Each result is dropped before the next pair. **Returned B** is the largest amount freed by dropping a path result, already included in peak; the two maxima can come from different pairs. On a cutoff miss, work can use heap without returning any. A zero peak means no heap allocation, not zero memory: stack, input buffers, allocator metadata and process RSS are not measured.

Both memory tables show each engine’s peak heap divided by Bitfront’s peak heap for the same row. `2.00×` means the engine used twice as much heap as Bitfront. A dash (`—`) means Bitfront used zero heap bytes, so no ratio is defined. Returned-path sizes stay in bytes.

The timing and memory tables show the unedited **Template bp**, the actual **Query bp** and **Target bp**, and the unbounded **Exact distance** (even for cutoff misses). Values that differ across pairs are shown as ranges. **Global extra edits** is the exact distance minus the absolute query–target length difference, computed separately for each pair before taking the range. It is shown as `—` for prefix and suffix searches.

On macOS, if `bindgen` cannot find Homebrew LLVM, set `LIBCLANG_PATH="$(brew --prefix llvm)/lib"`.

To time one workload, pass `--case MODE TEMPLATE_BP RATE_BP [CUTOFF]`. `MODE` is `global`, `prefix`, or `suffix`; `TEMPLATE_BP` is the positive unedited length in bases, and `RATE_BP` is the number of edits per 10,000 template bases (0–10000). Omit `CUTOFF` for an unbounded search. For example, this measures traceback on a 10,000 bp template at 1% edits against a 40,000 bp target:

```sh
sh tools/align-bench/run.sh --task path --case prefix 10000 100 --samples 50
```

Custom global cases can add trailing bases to the generated target with `--flank BP` (a nonnegative integer, default 0). Prefix and suffix cases always use a flank three times the template length, placed at the free end. To select the standard flanked global workload:

```sh
sh tools/align-bench/run.sh --task path --case global 256 500 --flank 768
```

Add `--swap` to a custom global case to exchange its query and target after generation, without changing the seeds. Reports mark these rows as `global (swapped)` and show lengths in the measured orientation. Swapped Criterion IDs add `-swap` before the optional `-path` suffix; unswapped IDs are unchanged. The `--flank` and `--swap` flags may precede or follow `--case`, but cannot be used with `--smoke` or with prefix or suffix cases.

Numeric bounded `--case` inputs replace selected template bases with `N`, which never occurs in the other input. The standard and smoke panels include hits at `CUTOFF=32` and misses at `31` for all three searches. These prefix and suffix cases use a 256 bp query and a 1,024 bp target.

The standard panel also includes a 10,000 bp global miss at `999` and hit at `1000` on the same four pairs. Run `sh tools/align-bench/run.sh --case global 10000 1000 999` (change the last argument for the hit) to get score and owned-path timings and their separate heap tables.

The standard panel (including CI, not `--smoke`) also runs three named 10,000 bp global controls for cutoff regressions. `mixed-below` and `mixed-at` use the **same four** shared-alphabet mixed-indel pairs. For each pair, `d` is its unbounded rust-edlib distance, computed before timing and heap sampling. Its cutoff is `d−1` for the miss or `d` for the hit: a cutoff range in the table is **not one shared bound**. `sparse-loose` has eight dispersed T→A substitutions on each pair and uses a fixed cutoff of `4096`. All three use the normal score, owned-path and heap reports. Run one without the rest of the panel:

```sh
sh tools/align-bench/run.sh --case mixed-below   # also: mixed-at, sparse-loose
```

## Compare with a previous Bitfront commit

Pass `--compare REV` first to add the specified revision to both the timing and memory tables for **score rows only**. Older revisions need not provide a path interface. `REV` is a commit, tag, or branch in this local Git repository. For example, compare the current working tree (including uncommitted changes) with the initial release:

```sh
RUSTUP_TOOLCHAIN=1.99.1 sh tools/align-bench/run.sh --compare 23f3d9c639a6a9950eabef29853dcdb2288c2180
```

Append `--smoke` or a custom `--case` (including `--flank` and `--swap` for global cases) to use a smaller panel. For example, `sh tools/align-bench/run.sh --compare HEAD --samples 50 --case suffix 256 1250 31` measures a bounded suffix miss. The `previous / Bitfront` column compares median time in the timing table and peak heap in the memory table; above 1× favors the working tree in both. Both versions score the same byte slices in the same release binary. The runner checks each result against the unbounded rust-edlib score for the same search mode before timing. Criterion warms up and measures one engine at a time.

This option requires Python 3 and Git. It resolves `REV` to a full commit ID and loads that source through a local Git dependency, without switching branches or modifying either tracked lockfile. Cargo resolves both versions together in a temporary project seeded with the benchmark lockfile; the chosen commit must have compatible dependencies and expose `distance`, `prefix_distance`, `suffix_distance`, and their `_at_most` variants.

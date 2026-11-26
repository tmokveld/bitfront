# Running the alignment benchmark

The benchmark compares Bitfront with direct calls to `rust-edlib` 0.3.1 from crates.io and with WFA2-lib. Both Rust scorers receive the same byte slices, including for suffix scoring. The direct rust-edlib baseline uses `align` with an inclusive `max_distance` for bounded global rows and its score-only functions for unbounded rows.

The [CI benchmark](../../.github/workflows/benchmark.yml) builds WFA2-lib with a C/C++ toolchain, Make, CMake 3.16 or newer, and `libclang`. Bitfront itself does not need these tools.

From the repository root, run the locked release benchmark used for the [README table](../../README.md#benchmark-results):

```sh
RUSTUP_TOOLCHAIN=1.99.1 sh tools/align-bench/run.sh
```

Each row uses generator indices `0..4`. Add `--smoke` for a smaller panel. On macOS, if `bindgen` cannot find Homebrew LLVM, set `LIBCLANG_PATH="$(brew --prefix llvm)/lib"`.

To time one workload, pass `--case MODE SOURCE_BP RATE_BP`. `MODE` is `global`, `prefix`, or `suffix`; `SOURCE_BP` is the source length in bases, and `RATE_BP` is the number of edits per 10,000 source bytes. For example, this scores a 10,000 bp source at 1% edits against a 40,000 bp target:

```sh
sh tools/align-bench/run.sh --case prefix 10000 100
```

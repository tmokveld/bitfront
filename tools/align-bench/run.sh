#!/bin/sh
set -eu

bench=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

if [ "${1:-}" = "--compare" ]; then
    shift
    case "${1:-}" in ''|-h|--help) exec python3 "$bench/compare.py" "$@" ;; esac
    set -- python3 "$bench/compare.py" "$@"
else
    set -- cargo bench --locked --manifest-path "$bench/Cargo.toml" --bench align -- "$@"
fi

# Keep stdout usable as Markdown (including CI's job summary). Criterion's
# native progress and confidence intervals stay visible on stderr.
CRITERION_HOME=${CRITERION_HOME:-${CARGO_TARGET_DIR:-$bench/target}/criterion}
case "$CRITERION_HOME" in /*) ;; *) CRITERION_HOME="$PWD/$CRITERION_HOME" ;; esac
export CRITERION_HOME
summary=$(mktemp "${TMPDIR:-/tmp}/bitfront-bench.XXXXXXXX")
trap 'rm -f "$summary"' 0
BITFRONT_BENCH_SUMMARY="$summary" "$@" >&2
cat "$summary"

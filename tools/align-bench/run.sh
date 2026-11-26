#!/bin/sh
set -eu

exec cargo run --release --locked --manifest-path "$(dirname "$0")/Cargo.toml" -- "$@"

#!/usr/bin/env bash
set -euo pipefail

echo "========================================"
echo "  Float vs CT Physics Sim — Benchmarks"
echo "========================================"
echo

# Release build
echo "Building release binary…"
cargo build --release 2>&1

echo
echo "--- Release run ---"
time cargo run --release --quiet

echo
echo "--- Debug run (slower, for comparison) ---"
echo "Building debug binary…"
cargo build 2>&1
time cargo run --quiet

echo
echo "Done."

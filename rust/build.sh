#! /usr/bin/bash
rm -rf target
cargo build --release --features flat
cp target/release/gils_rvnd solve

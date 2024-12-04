#! /usr/bin/bash
rm -rf target
cargo build --release --features flat -j8
cp target/release/gils_rvnd solve

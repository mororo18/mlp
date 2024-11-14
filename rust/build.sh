#! /usr/bin/bash
cargo build --release --features flat
cp target/debug/main solve

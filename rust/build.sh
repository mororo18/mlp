#! /usr/bin/bash
rm main
cargo build --release
cp target/release/main .

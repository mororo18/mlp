#! /usr/bin/bash
rm main
cargo build
cp target/debug/main .

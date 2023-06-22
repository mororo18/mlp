#! /bin/bash

rm target/debug/main
cargo build
cp target/debug/main .

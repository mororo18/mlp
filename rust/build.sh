#! /usr/bin/bash
rustc -L dependencies/rand-0.8.4/target/debug/deps/ main.rs && ./main.rs
rm main

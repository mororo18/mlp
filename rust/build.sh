#! /usr/bin/bash
rustc -C opt-level=3 -g -L dependencies/rand-0.8.4/target/debug/deps/ main.rs && ./main
#rm main

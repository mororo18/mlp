#!/usr/bin/env bash
# Diagnóstico do P7 (lin318 lento no JavaScript): o nosso JavaScript no gr202 e
# no lin318 com Node 20, 22 e 24, no servidor, com as mesmas matrizes de
# distância do bench.py (.bench/dm_cache).
#
# Uso (no servidor): NIX_IMAGE_DIGEST=... llm-translations/diag/node_versions.sh <saida>
# Node 20 e 22 vêm do nixos-25.05; o 24 é o do flake.nix das rodadas.
set -euo pipefail

LT=$(cd "$(dirname "$0")/.." && pwd)
REPO=$(dirname "$LT")
OUT=$1
IMAGE=docker.io/nixos/nix@sha256:${NIX_IMAGE_DIGEST:?set NIX_IMAGE_DIGEST}
NIX="nix --extra-experimental-features nix-command --extra-experimental-features flakes"

if [ -z "${INSIDE:-}" ]; then
    mkdir -p "$OUT/flake"
    OUT=$(cd "$OUT" && pwd)
    git -C "$REPO" rev-parse HEAD > "$OUT/commit.txt"
    cp "$REPO/flake.nix" "$REPO/flake.lock" "$OUT/flake/"
    taskset -c 8,9,10,11 podman run --rm -v fa-nix:/nix -v "$REPO":/repo -v "$OUT":/out -w /repo \
        -e INSIDE=1 -e NIX_IMAGE_DIGEST="$NIX_IMAGE_DIGEST" "$IMAGE" \
        bash /repo/llm-translations/diag/node_versions.sh /out
    echo FIM
    exit 0
fi

run_pair() {
    local tag=$1; shift
    for inst in gr202 lin318; do
        local box=/tmp/box_${tag}_${inst}
        rm -rf "$box"; mkdir -p "$box"
        cp -r /repo/javascript "$box/js"
        cp "/repo/.bench/dm_cache/${inst}.tsp.distance_matrix" "$box/distance_matrix"
        ( cd "$box/js" && "$@" bash -c 'node --version; node main.js' > "/out/node${tag}_${inst}.txt" 2>&1 ) &
    done
    wait
}

run_pair 20 $NIX shell github:NixOS/nixpkgs/nixos-25.05#nodejs_20 -c
run_pair 22 $NIX shell github:NixOS/nixpkgs/nixos-25.05#nodejs_22 -c
run_pair 24 $NIX develop path:/out/flake -c

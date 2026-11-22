#!/usr/bin/env bash
# Diagnóstico do P16 (dois patamares do Node 24): o nosso JavaScript numa
# instância, várias execuções em paralelo, com conjuntos de flags do V8
# intercalados, no servidor, com o Node 24 do flake.nix.
#
# Uso (no servidor): NIX_IMAGE_DIGEST=... [INSTANCE=pr299] [REPS=10] [JOBS=4] \
#     [FLAGSETS="default|--single-threaded|--no-concurrent-recompilation"] \
#     llm-translations/diag/node_flags.sh <saida>
# Cada execução grava <saida>/<n>_<rep>.txt (n = posição do conjunto em
# FLAGSETS) com as flags, o tempo de relógio e a saída do programa.
set -euo pipefail

LT=$(cd "$(dirname "$0")/.." && pwd)
REPO=$(dirname "$LT")
OUT=$1
IMAGE=docker.io/nixos/nix@sha256:${NIX_IMAGE_DIGEST:?set NIX_IMAGE_DIGEST}
NIX="nix --extra-experimental-features nix-command --extra-experimental-features flakes"
DEFAULT_SETS="default|--single-threaded|--no-concurrent-recompilation"

if [ -z "${INSIDE:-}" ]; then
    mkdir -p "$OUT/flake"
    OUT=$(cd "$OUT" && pwd)
    git -C "$REPO" rev-parse HEAD > "$OUT/commit.txt"
    cp "$REPO/flake.nix" "$REPO/flake.lock" "$OUT/flake/"
    taskset -c 8,9,10,11 podman run --rm -v fa-nix:/nix -v "$REPO":/repo -v "$OUT":/out -w /repo \
        -e INSIDE=1 -e NIX_IMAGE_DIGEST="$NIX_IMAGE_DIGEST" -e INSTANCE="${INSTANCE:-pr299}" \
        -e REPS="${REPS:-10}" -e JOBS="${JOBS:-4}" -e FLAGSETS="${FLAGSETS:-$DEFAULT_SETS}" "$IMAGE" \
        bash /repo/llm-translations/diag/node_flags.sh /out
    echo FIM
    exit 0
fi

DM=/repo/.bench/dm_cache/$INSTANCE.tsp.distance_matrix
if [ ! -f "$DM" ]; then
    $NIX develop path:/out/flake -c bash -c "cd /repo/mlp-instances && ./load $INSTANCE.tsp > /dev/null && cp ../distance_matrix $DM"
fi

IFS='|' read -ra SETS <<< "$FLAGSETS"
$NIX develop path:/out/flake -c bash -c 'node --version; for f in "$@"; do [ "$f" = default ] || node $f -e 0 && echo "ok $f"; done' _ "${SETS[@]}" > /out/flags.txt 2>&1

run_one() {
    local n=$1 rep=$2 flags=$3 box=/tmp/box_${1}_${2}
    [ "$flags" = default ] && flags=
    rm -rf "$box"; mkdir -p "$box"
    cp -r /repo/javascript "$box/js"
    cp "$DM" "$box/distance_matrix"
    ( cd "$box/js" && $NIX develop path:/out/flake -c bash -c \
        'echo "flags: $1"; s=$(date +%s%N); node $1 main.js; e=$(date +%s%N); echo "wall_ms: $(( (e - s) / 1000000 ))"' _ "$flags" \
        > "/out/${n}_${rep}.txt" 2>&1 )
    rm -rf "$box"
}

for rep in $(seq 0 $((REPS - 1))); do
    for n in "${!SETS[@]}"; do
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
        run_one "$n" "$rep" "${SETS[$n]}" &
    done
done
wait

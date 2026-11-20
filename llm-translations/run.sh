#!/usr/bin/env bash
# Rodada do experimento no servidor (q-dual).
#
# Uso: llm-translations/run.sh <saida> "<linguagens nossas>" <variante> [<variante> ...]
#   ex.: llm-translations/run.sh ~/mlp-test/fa/runs/base "java julia cplusplus" cpp_java.base rust_cpp.base
#
# Fora do container: confere que o código das nossas implementações e desta
# pasta não tem mudanças, grava o commit, liga o log de carga (uma leitura por
# minuto: hora, carga de 1/5/15 min, processos acima de 50% de CPU) e abre o
# container com o ambiente do flake.nix deste repositório, preso aos CPUs 8 a
# 11 (taskset; o bench.py fixa cada execução num desses núcleos).
# Dentro do container: grava as versões das toolchains, monta as variantes e
# roda o bench.py com -n $REPS -j 4 nas instâncias de $INST (padrão: 10 e
# instances.txt; o teste de build usa REPS=1 INST=smoke.txt).
set -euo pipefail

LT=$(cd "$(dirname "$0")" && pwd)
REPO=$(dirname "$LT")
OUT=$1
LANGS=$2
shift 2
IDS="$*"
REPS=${REPS:-10}
INST=${INST:-instances.txt}
IMAGE=docker.io/nixos/nix@sha256:${NIX_IMAGE_DIGEST:?set NIX_IMAGE_DIGEST}

if [ -z "${IN_NIX_SHELL:-}" ]; then
    LANG_DIRS="c cplusplus csharp fortran go java javascript julia lua python rust swift"
    if ! git -C "$REPO" diff --quiet HEAD -- $LANG_DIRS llm-translations; then
        echo "uncommitted changes in the implementations or in llm-translations" >&2
        exit 1
    fi
    mkdir -p "$OUT/flake"
    OUT=$(cd "$OUT" && pwd)
    git -C "$REPO" rev-parse HEAD > "$OUT/commit.txt"
    echo "$IMAGE" > "$OUT/image.txt"
    cp "$REPO/flake.nix" "$REPO/flake.lock" "$OUT/flake/"
    ( while true; do
        echo "$(date +%H:%M) $(cut -d' ' -f1-3 /proc/loadavg) $(ps -eo pcpu= | awk '$1>50' | wc -l)"
        sleep 60
      done ) >> "$OUT/carga.log" &
    LOGGER=$!
    trap 'kill $LOGGER 2>/dev/null || true' EXIT
    taskset -c 8,9,10,11 podman run --rm \
        -v fa-nix:/nix -v "$REPO":/repo -v "$OUT":/out -w /repo -e REPS="$REPS" -e INST="$INST" \
        -e NIX_IMAGE_DIGEST="${NIX_IMAGE_DIGEST}" \
        "$IMAGE" nix --extra-experimental-features "nix-command flakes" \
        develop path:/out/flake -c bash /repo/llm-translations/run.sh /out "$LANGS" $IDS
    echo FIM
    exit 0
fi

{
    g++ --version | head -1
    gfortran --version | head -1
    go version
    javac -version 2>&1
    java -version 2>&1 | head -1
    node --version
    swift --version 2>&1 | head -1
    rustc --version
    julia --version
    python3 --version
    echo "affinity: $(python3 -c 'import os; print(sorted(os.sched_getaffinity(0)))')"
} > "$OUT/versions.txt"

mkdir -p "$REPO/benchmarks"
cp "$LT/harness/bench.py" "$LT/harness/core.py" "$REPO/benchmarks/"
python3 "$LT/make_variants.py" --out "$OUT/work" --ids $IDS
python3 "$REPO/benchmarks/bench.py" -l $LANGS -t $IDS --translations-dir "$OUT/work" \
    -n "$REPS" -j 4 -I "$LT/$INST" --out "$OUT/csv"

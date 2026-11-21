"""Gaps de uma rodada do experimento, a partir dos CSVs brutos do bench.py.

Para cada variante (coluna `source` = id do manifesto), compara com a nossa
implementação na linguagem de destino da tradução (C++ para as traduções
`<ling>_cpp`): gap = 100 x (média, nas instâncias, da razão variante / nossa
- 1), com a razão de cada instância calculada sobre as médias das execuções.
Também dá o intervalo das razões por instância e confere o custo: toda
execução de toda variante tem de reproduzir o custo da nossa implementação.

Uso: python analyze.py <pasta csv da rodada> [--sem berlin52.tsp]
"""
import argparse
import csv
import glob
import os
import sys
from statistics import mean, pstdev

OURS = {"cpp_java": "java", "cpp_julia": "julia", "cpp_swift": "swift", "cpp_go": "go",
        "cpp_fortran": "fortran", "cpp_javascript": "javascript"}


def our_lang(variant):
    pair = variant.split(".")[0]
    if pair.startswith("ours_"):
        return pair[len("ours_"):]
    return OURS.get(pair, "cplusplus")


def load(folder):
    runs = {}
    for path in glob.glob(os.path.join(folder, "*.csv")):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                runs.setdefault(r["source"], {}).setdefault(r["instance"], []).append(r)
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--sem", nargs="*", default=[])
    args = ap.parse_args()
    runs = load(args.folder)
    ok = True
    print(f"{'variant':24} {'ours':11} {'n':>3} {'gap %':>7} {'min':>6} {'max':>6} {'max cv %':>8}  cost")
    for v in sorted(k for k in runs if "." in k):
        ref = runs.get(our_lang(v))
        if not ref:
            print(f"{v:24} missing our implementation {our_lang(v)}")
            ok = False
            continue
        ratios, cvs, cost_ok = [], [], True
        for inst, rs in sorted(runs[v].items()):
            if inst in args.sem or inst not in ref:
                continue
            ref_costs = {float(r["COST"]) for r in ref[inst] if r["timed_out"] != "True"}
            costs = {float(r["COST"]) for r in rs if r["timed_out"] != "True"}
            if len(ref_costs) != 1 or costs != ref_costs or any(r["timed_out"] == "True" for r in rs):
                cost_ok = False
            t = [float(r["TIME"]) for r in rs]
            tr = [float(r["TIME"]) for r in ref[inst]]
            ratios.append(mean(t) / mean(tr))
            cvs.append(100 * pstdev(t) / mean(t))
        ok &= cost_ok
        print(f"{v:24} {our_lang(v):11} {len(ratios):3d} {100 * (mean(ratios) - 1):+7.1f} "
              f"{min(ratios):6.3f} {max(ratios):6.3f} {max(cvs):8.1f}  {'ok' if cost_ok else 'MISMATCH'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

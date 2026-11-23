"""Comparações da análise de fatores (Seções 7.5 e 7.6 do artigo).

Cada comparação é a razão do tempo médio de dois programas medidos na mesma
rodada (runs/<rodada>/csv/<programa>.csv), por instância. Gap = 100 x (média
das razões nas instâncias - 1), o método da Figura 7; também o gap sem o
berlin52. Colunas das tabelas:
- gap: tradução / nossa implementação;
- ours_in_translation: tradução com a nossa escolha / nossa implementação;
- translation_in_ours: tradução / nossa implementação com a escolha da tradução;
- both: tradução / nossa implementação, com a mesma mudança nos dois
  (aquecimento do JIT);
- control: variante da nossa implementação / nossa implementação.

Uso: python tables.py [--out comparisons.csv]
"""
import argparse
import csv
from pathlib import Path
from statistics import mean, stdev

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
INSTANCES = ["berlin52", "kroA100", "gr202", "lin318"]

C = [
    ("pitfalls", "to", "Rust -> C++", "", "gap", "factors-20261121", "rust_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Rust -> C++", "Distance matrix", "ours_in_translation", "factors-20261121", "rust_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "C# -> C++", "", "gap", "factors-20261121", "csharp_cpp.base", "cplusplus"),
    ("pitfalls", "to", "C# -> C++", "Distance matrix", "ours_in_translation", "factors-20261121", "csharp_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Swift -> C++", "", "gap", "factors-20261121", "swift_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Swift -> C++", "Distance matrix", "ours_in_translation", "factors-20261121", "swift_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Swift -> C++", "Pass by value", "ours_in_translation", "factors-20261121", "swift_cpp.route_copy", "cplusplus"),
    ("pitfalls", "to", "Go -> C++", "", "gap", "factors-20261121", "go_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Go -> C++", "Distance matrix", "ours_in_translation", "factors-20261121", "go_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Fortran -> C++", "", "gap", "factors-20261121", "fortran_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Fortran -> C++", "Distance matrix", "ours_in_translation", "factors-20261121", "fortran_cpp.row_matrix", "cplusplus"),

    ("pitfalls", "from", "C++ -> JavaScript", "", "gap", "jsf64-20261122", "cpp_javascript.base", "javascript"),
    ("pitfalls", "from", "C++ -> JavaScript", "Subsequence storage", "ours_in_translation", "jsf64-20261122", "cpp_javascript.array", "javascript"),
    ("pitfalls", "from", "C++ -> JavaScript", "Subsequence storage", "translation_in_ours", "jsf64-20261122", "cpp_javascript.base", "ours_javascript.float64array"),
    ("pitfalls", "from", "C++ -> Swift", "", "gap", "factors-20261121", "cpp_swift.base", "swift"),
    ("pitfalls", "from", "C++ -> Swift", "Copy of the solution", "translation_in_ours", "factors-20261121", "cpp_swift.base", "ours_swift.route_copy"),
    ("pitfalls", "from", "C++ -> Swift", "Copy of the solution", "ours_in_translation", "p3p8-20261121", "cpp_swift.full_copy", "swift"),
    ("pitfalls", "from", "C++ -> Java", "", "gap", "p11p12-20261121", "cpp_java.base", "java"),
    ("pitfalls", "from", "C++ -> Java", "Subsequence storage", "ours_in_translation", "java-factors-20261121", "cpp_java.aos", "java"),
    ("pitfalls", "from", "C++ -> Java", "JIT warm-up", "both", "p11p12-20261121", "cpp_java.warm", "ours_java.warm"),
    ("pitfalls", "from", "C++ -> Julia", "", "gap", "p11p12-20261121", "cpp_julia.base", "julia"),
    ("pitfalls", "from", "C++ -> Julia", "JIT compilation", "both", "p11p12-20261121", "cpp_julia.warm", "ours_julia.warm"),

    ("other", "to", "Go -> C++", "Epsilon", "ours_in_translation", "factors-20261121", "go_cpp.epsilon", "cplusplus"),
    ("other", "to", "Go -> C++", "Epsilon + distance matrix", "ours_in_translation", "combo-20261121", "go_cpp.epsilon_row_matrix", "cplusplus"),
    ("other", "to", "Fortran -> C++", "Epsilon", "ours_in_translation", "factors-20261121", "fortran_cpp.epsilon", "cplusplus"),
    ("other", "to", "Fortran -> C++", "Epsilon + distance matrix", "ours_in_translation", "combo-20261121", "fortran_cpp.epsilon_row_matrix", "cplusplus"),
    ("other", "to", "Rust -> C++", "Index type", "ours_in_translation", "factors-20261121", "rust_cpp.index32", "cplusplus"),
    ("other", "to", "C# -> C++", "Index type", "ours_in_translation", "factors-20261121", "csharp_cpp.index32", "cplusplus"),
    ("other", "to", "Swift -> C++", "Index type", "ours_in_translation", "factors-20261121", "swift_cpp.index32", "cplusplus"),
    ("other", "to", "Go -> C++", "Index type", "ours_in_translation", "factors-20261121", "go_cpp.index32", "cplusplus"),
    ("other", "to", "Fortran -> C++", "Index type", "ours_in_translation", "factors-20261121", "fortran_cpp.index32", "cplusplus"),
    ("other", "to", "Rust -> C++", "Index type", "translation_in_ours", "factors-20261121", "rust_cpp.base", "ours_cplusplus.index64"),
    ("other", "to", "C# -> C++", "Index type", "translation_in_ours", "factors-20261121", "csharp_cpp.base", "ours_cplusplus.index64"),
    ("other", "to", "Swift -> C++", "Index type", "translation_in_ours", "factors-20261121", "swift_cpp.base", "ours_cplusplus.index64"),
    ("other", "to", "Go -> C++", "Index type", "translation_in_ours", "factors-20261121", "go_cpp.base", "ours_cplusplus.index64"),
    ("other", "to", "Fortran -> C++", "Index type", "translation_in_ours", "factors-20261121", "fortran_cpp.base", "ours_cplusplus.index64"),
    ("other", "from", "C++ -> Java", "Accessors + local copies", "translation_in_ours", "p11p12-20261121", "cpp_java.base", "ours_java.direct_no_locals"),
    ("other", "from", "C++ -> Java", "Accessors + local copies + subsequence storage", "translation_in_ours", "p11p12-20261121", "cpp_java.base", "ours_java.soa"),
    ("other", "from", "C++ -> Java", "All, with JIT warm-up", "translation_in_ours", "p11p12-20261121", "cpp_java.warm", "ours_java.soa_warm"),
    ("other", "from", "C++ -> Julia", "1-based city index, with JIT warm-up", "translation_in_ours", "p11p12-20261121", "cpp_julia.warm", "ours_julia.zero_based_warm"),
    ("other", "from", "C++ -> Go", "", "gap", "factors-20261121", "cpp_go.base", "go"),
    ("other", "from", "C++ -> Go", "Accessor", "translation_in_ours", "factors-20261121", "cpp_go.base", "ours_go.direct_index"),
    ("other", "from", "C++ -> Go", "Accessor", "ours_in_translation", "p3p8-20261121", "cpp_go.at_accessor", "go"),

    ("controls", "from", "JavaScript", "Float64Array instead of Array", "control", "jsf64-20261122", "ours_javascript.float64array", "javascript"),
    ("controls", "from", "C++ -> JavaScript", "Array instead of Float64Array, against the translation", "control", "jsf64-20261122", "cpp_javascript.array", "cpp_javascript.base"),
    ("controls", "from", "C++ -> JavaScript", "Both Float64Array: translation / ours", "control", "p20-20261123", "cpp_javascript.base", "javascript"),
    ("controls", "from", "C++ -> JavaScript", "Both Float64Array: interleaved translation / ours", "control", "p20-20261123", "cpp_javascript.interleaved", "javascript"),
    ("controls", "from", "C++ -> JavaScript", "Both Float64Array: translation / ours with three arrays", "control", "p20-20261123", "cpp_javascript.base", "ours_javascript.three_arrays"),
    ("controls", "from", "JavaScript", "Float64Array: three arrays instead of interleaved", "control", "p20-20261123", "ours_javascript.three_arrays", "javascript"),
    ("controls", "from", "C++ -> JavaScript", "Float64Array: interleaved instead of three arrays, against the translation", "control", "p20-20261123", "cpp_javascript.interleaved", "cpp_javascript.base"),
    ("controls", "to", "C++", "size_t counter truncated to 32 bits", "control", "p3p8-20261121", "ours_cplusplus.size_t_counter", "cplusplus"),
    ("controls", "to", "C++", "64-bit index", "control", "factors-20261121", "ours_cplusplus.index64", "cplusplus"),
    ("controls", "to", "C# -> C++", "", "gap", "p14-20261121", "csharp_cpp.base", "cplusplus"),
    ("controls", "to", "C# -> C++", "Index type, subsequences only", "ours_in_translation", "p14-20261121", "csharp_cpp.index32_seq", "cplusplus"),
    ("controls", "to", "C# -> C++", "Index type, distance matrix only", "ours_in_translation", "p14-20261121", "csharp_cpp.index32_matrix", "cplusplus"),
    ("controls", "to", "C# -> C++", "Index type, subsequences, non-negative", "ours_in_translation", "p14-20261121", "csharp_cpp.index32_seq_assume", "cplusplus"),
    ("controls", "to", "Rust -> C++", "Index type computed in int", "ours_in_translation", "p3p8-20261121", "rust_cpp.index32_int", "cplusplus"),
]


def reference_costs():
    costs = {}
    with open(RUNS / "base-20261120" / "csv" / "cplusplus.csv", newline="") as f:
        for r in csv.DictReader(f):
            costs[r["instance"].removesuffix(".tsp")] = r["COST"]
    return costs


def load(round_, prog, costs):
    times = {}
    with open(RUNS / round_ / "csv" / f"{prog}.csv", newline="") as f:
        for r in csv.DictReader(f):
            inst = r["instance"].removesuffix(".tsp")
            assert r["COST"] == costs[inst], (round_, prog, inst, r["COST"])
            times.setdefault(inst, []).append(float(r["TIME"]))
    return times


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "comparisons.csv"))
    args = ap.parse_args()
    costs = reference_costs()
    rows = []
    for table, direction, pair, factor, column, round_, num, den in C:
        a, b = load(round_, num, costs), load(round_, den, costs)
        ratios = [mean(a[i]) / mean(b[i]) for i in INSTANCES]
        cv = max(stdev(v) / mean(v) for t in (a, b) for v in (t[i] for i in INSTANCES))
        rows.append({
            "table": table, "direction": direction, "pair": pair, "factor": factor, "column": column,
            "round": round_, "numerator": num, "denominator": den,
            **{f"ratio_{i}": f"{r:.4f}" for i, r in zip(INSTANCES, ratios)},
            "gap_pct": f"{100 * (mean(ratios) - 1):.1f}",
            "gap_no_berlin52_pct": f"{100 * (mean(ratios[1:]) - 1):.1f}",
            "max_cv_pct": f"{100 * cv:.1f}",
        })
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['table']:9} {r['pair']:18} {r['column']:20} {r['gap_pct']:>6} {r['gap_no_berlin52_pct']:>6}  {r['factor']}")


if __name__ == "__main__":
    main()

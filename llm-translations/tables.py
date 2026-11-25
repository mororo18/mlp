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
# C++ -> Lua e C++ -> Python no Lua 5.2 e no CPython (rodadas *-interp-*).
INTERP = ["berlin52", "kroA100", "lin105"]
ALL = INSTANCES + ["lin105"]
B4C, B4I, B5C = "b4-comp-20261124", "b4-interp-20261124", "b5-comp-20261125"
B7C, B7S, B7I = "b7-comp-20261125", "b7-short-20261125", "b7-interp-20261125"
B8I, B9C, B9I = "b8-lua-20261125", "b9-comp-20261125", "b9-interp-20261125"
SHORT = ["berlin52", "kroA100"]


def instances(round_):
    if "short" in round_:
        return SHORT
    return INTERP if "interp" in round_ or "-lua-" in round_ else INSTANCES

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

    # Frente B (PLANO_FECHAMENTO.md): pares que só tinham a análise antiga.
    ("pitfalls", "from", "C++ -> Rust", "", "gap", B4C, "cpp_rust.base", "rust"),
    ("pitfalls", "from", "C++ -> Rust", "Bounds checks", "ours_in_translation", B4C, "cpp_rust.unchecked", "rust"),
    ("pitfalls", "from", "C++ -> Rust", "Distance matrix", "ours_in_translation", B4C, "cpp_rust.flat_matrix", "rust"),
    ("pitfalls", "from", "C++ -> Rust", "Bounds checks + distance matrix", "ours_in_translation", B4C, "cpp_rust.unchecked_flat_matrix", "rust"),
    ("pitfalls", "from", "C++ -> C#", "", "gap", B4C, "cpp_csharp.base", "csharp"),
    ("pitfalls", "from", "C++ -> C#", "Bounds checks", "ours_in_translation", B4C, "cpp_csharp.unchecked_seq", "csharp"),
    ("pitfalls", "from", "C++ -> C#", "Distance matrix", "ours_in_translation", B4C, "cpp_csharp.flat_matrix", "csharp"),
    ("pitfalls", "from", "C++ -> C#", "Bounds checks + distance matrix", "ours_in_translation", B4C, "cpp_csharp.unchecked_seq_flat_matrix", "csharp"),
    ("pitfalls", "from", "C++ -> LuaJIT", "", "gap", B4C, "cpp_lua.luajit", "lua"),
    ("pitfalls", "from", "C++ -> LuaJIT", "Subsequence storage", "ours_in_translation", B4C, "cpp_lua.luajit_subseq_array", "lua"),
    ("pitfalls", "from", "C++ -> PyPy", "", "gap", B4C, "cpp_pypy.base", "python"),
    ("pitfalls", "from", "C++ -> PyPy", "Element types", "ours_in_translation", B4C, "cpp_pypy.float_w", "python"),
    ("pitfalls", "from", "C++ -> PyPy", "Subsequence storage", "ours_in_translation", B4C, "cpp_pypy.interleaved", "python"),
    ("pitfalls", "from", "C++ -> PyPy", "Element types + subsequence storage", "ours_in_translation", B4C, "cpp_pypy.float_w_interleaved", "python"),
    ("pitfalls", "from", "C++ -> Lua", "", "gap", B4I, "cpp_lua.lua", "lua"),
    ("pitfalls", "from", "C++ -> Lua", "Subsequence storage", "ours_in_translation", B4I, "cpp_lua.lua_subseq_array", "lua"),
    ("pitfalls", "from", "C++ -> Python", "", "gap", B4I, "cpp_python.base", "python"),
    ("pitfalls", "from", "C++ -> Python", "Subsequence storage", "ours_in_translation", B4I, "cpp_python.nested_lists", "python"),
    ("pitfalls", "from", "C++ -> Python", "Element types", "ours_in_translation", B4I, "cpp_python.float_w", "python"),
    ("pitfalls", "to", "C -> C++", "", "gap", B5C, "c_cpp.base", "cplusplus"),
    ("pitfalls", "to", "C -> C++", "Distance matrix", "ours_in_translation", B5C, "c_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Lua -> C++", "", "gap", B5C, "lua_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Lua -> C++", "Distance matrix", "ours_in_translation", B5C, "lua_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Python -> C++", "", "gap", B4C, "python_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Python -> C++", "Subsequence storage", "ours_in_translation", B4C, "python_cpp.flat_seq", "cplusplus"),
    ("pitfalls", "to", "Python -> C++", "Distance matrix", "ours_in_translation", B4C, "python_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Java -> C++", "", "gap", B4C, "java_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Java -> C++", "Subsequence storage", "ours_in_translation", B4C, "java_cpp.flat_seq", "cplusplus"),
    ("pitfalls", "to", "Java -> C++", "Distance matrix", "ours_in_translation", B4C, "java_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "JavaScript -> C++", "", "gap", B4C, "javascript_cpp.base", "cplusplus"),
    ("pitfalls", "to", "JavaScript -> C++", "Distance matrix", "ours_in_translation", B4C, "javascript_cpp.row_matrix", "cplusplus"),
    ("pitfalls", "to", "Julia -> C++", "", "gap", B5C, "julia_cpp.base", "cplusplus"),
    ("pitfalls", "to", "Julia -> C++", "Distance matrix", "ours_in_translation", B5C, "julia_cpp.row_matrix", "cplusplus"),

    ("other", "from", "C++ -> C", "", "gap", "b-comp-20261124", "cpp_c.base", "c"),
    ("other", "from", "C++ -> LuaJIT", "Function calls", "ours_in_translation", B4C, "cpp_lua.luajit_inline_to_1d", "lua"),
    ("other", "from", "C++ -> LuaJIT", "Fields re-read in loops", "ours_in_translation", B4C, "cpp_lua.luajit_local_fields", "lua"),
    ("other", "from", "C++ -> LuaJIT", "All", "ours_in_translation", B4C, "cpp_lua.luajit_all", "lua"),
    ("other", "from", "C++ -> Lua", "Function calls", "ours_in_translation", B4I, "cpp_lua.lua_inline_to_1d", "lua"),
    ("other", "from", "C++ -> Lua", "Fields re-read in loops", "ours_in_translation", B4I, "cpp_lua.lua_local_fields", "lua"),
    ("other", "from", "C++ -> Lua", "Function calls + fields + subsequence storage", "ours_in_translation", B4I, "cpp_lua.lua_all", "lua"),
    ("other", "from", "C++ -> Python", "Function calls", "ours_in_translation", B4I, "cpp_python.inline_to_1d", "python"),
    ("other", "from", "C++ -> Python", "Fields re-read in loops", "ours_in_translation", B4I, "cpp_python.local_fields", "python"),
    ("other", "from", "C++ -> Python", "All", "ours_in_translation", B4I, "cpp_python.all", "python"),
    ("other", "to", "Lua -> C++", "Epsilon", "ours_in_translation", B5C, "lua_cpp.epsilon", "cplusplus"),
    ("other", "to", "Lua -> C++", "Epsilon + distance matrix", "ours_in_translation", B5C, "lua_cpp.epsilon_row_matrix", "cplusplus"),
    ("other", "to", "Python -> C++", "Epsilon", "ours_in_translation", B4C, "python_cpp.epsilon", "cplusplus"),
    ("other", "to", "Python -> C++", "Epsilon + subsequence storage + distance matrix", "ours_in_translation", B4C, "python_cpp.all", "cplusplus"),
    ("other", "to", "Java -> C++", "Epsilon", "ours_in_translation", B4C, "java_cpp.epsilon", "cplusplus"),
    ("other", "to", "Java -> C++", "Epsilon + subsequence storage + distance matrix", "ours_in_translation", B5C, "java_cpp.all", "cplusplus"),
    ("other", "to", "JavaScript -> C++", "Offsets", "ours_in_translation", B4C, "javascript_cpp.const_offsets", "cplusplus"),
    ("other", "to", "JavaScript -> C++", "Offsets + distance matrix", "ours_in_translation", B4C, "javascript_cpp.const_offsets_row_matrix", "cplusplus"),
    ("other", "to", "Julia -> C++", "Offsets", "ours_in_translation", B5C, "julia_cpp.const_offsets", "cplusplus"),
    ("other", "to", "Julia -> C++", "Offsets + distance matrix", "ours_in_translation", B5C, "julia_cpp.const_offsets_row_matrix", "cplusplus"),

    # Tabela 13 (tab:translation-closing): a combinação que fecha cada par.
    ("closing", "from", "C++ -> JavaScript", "", "gap", B5C, "cpp_javascript.base", "javascript"),
    ("closing", "from", "C++ -> JavaScript", "Array, interleaved, packed matrix (translation); constant positions, local copies (ours)", "both", B5C, "cpp_javascript.interleaved_array_packed", "ours_javascript.local_info_consts"),
    ("closing", "from", "C++ -> JavaScript", "Array, interleaved (translation); constant positions, local copies (ours)", "both", B5C, "cpp_javascript.interleaved_array", "ours_javascript.local_info_consts"),
    ("closing", "from", "C++ -> JavaScript", "Packed matrix", "ours_in_translation", B5C, "cpp_javascript.packed_matrix", "javascript"),
    ("closing", "from", "C++ -> JavaScript", "Array, interleaved", "ours_in_translation", B5C, "cpp_javascript.interleaved_array", "javascript"),
    ("closing", "from", "C++ -> JavaScript", "Array, interleaved, packed matrix", "ours_in_translation", B5C, "cpp_javascript.interleaved_array_packed", "javascript"),
    ("closing", "from", "JavaScript", "Constant positions, local copies", "control", B5C, "ours_javascript.local_info_consts", "javascript"),
    # B7: escolhas das traduções nas nossas (Java, LuaJIT, Python, PyPy) e aquecimento (C#, PyPy).
    ("closing", "from", "C++ -> Java", "", "gap", B7C, "cpp_java.base", "java"),
    ("closing", "from", "C++ -> Java", "Accessors + local copies + subsequence storage", "translation_in_ours", B7C, "cpp_java.base", "ours_java.soa"),
    ("closing", "from", "C++ -> Java", "Reading of the instance", "translation_in_ours", B7C, "cpp_java.base", "ours_java.bufread"),
    ("closing", "from", "C++ -> Java", "Accessors + local copies + subsequence storage + reading of the instance", "translation_in_ours", B7C, "cpp_java.base", "ours_java.soa_bufread"),
    ("closing", "from", "C++ -> LuaJIT", "Subsequence storage", "ours_in_translation", B7C, "cpp_lua.luajit_subseq_array", "lua"),
    ("closing", "from", "C++ -> LuaJIT", "Subsequence storage (translation); literal positions (ours)", "both", B7C, "cpp_lua.luajit_subseq_array", "ours_lua.literal_offsets"),
    ("closing", "from", "C++ -> LuaJIT", "Subsequence storage (translation); literal positions, loop locals (ours)", "both", B7C, "cpp_lua.luajit_subseq_array", "ours_lua.literal_offsets_loop_locals"),
    ("closing", "from", "LuaJIT", "Literal positions", "control", B7C, "ours_lua.literal_offsets", "lua"),
    ("closing", "from", "C++ -> PyPy", "Element types + subsequence storage", "ours_in_translation", B7C, "cpp_pypy.float_w_interleaved", "python"),
    ("closing", "from", "C++ -> PyPy", "Element types + subsequence storage (translation); inline index (ours)", "both", B7C, "cpp_pypy.float_w_interleaved", "ours_python.pypy_inline"),
    ("closing", "from", "PyPy", "Inline index", "control", B7C, "ours_python.pypy_inline", "python"),
    ("closing", "from", "C++ -> Python", "All (translation)", "ours_in_translation", B7I, "cpp_python.all", "python"),
    ("closing", "from", "C++ -> Python", "All (translation); local copies in the update (ours)", "both", B7I, "cpp_python.all", "ours_python.update_locals"),
    ("closing", "from", "C++ -> Python", "All (translation); local copies in the update, literal positions (ours)", "both", B7I, "cpp_python.all", "ours_python.update_locals_literals"),
    ("closing", "from", "C++ -> C#", "Bounds checks + distance matrix", "ours_in_translation", B7S, "cpp_csharp.unchecked_seq_flat_matrix", "csharp"),
    ("closing", "from", "C++ -> C#", "Bounds checks + distance matrix, JIT warm-up in both", "both", B7S, "cpp_csharp.unchecked_seq_flat_matrix_warm", "ours_csharp.warm"),
    ("closing", "from", "C#", "JIT warm-up", "control", B7S, "ours_csharp.warm", "csharp"),
    ("closing", "from", "C++ -> PyPy", "Element types + subsequence storage, short round", "ours_in_translation", B7S, "cpp_pypy.float_w_interleaved", "python"),
    ("closing", "from", "C++ -> PyPy", "Element types + subsequence storage, JIT warm-up in both", "both", B7S, "cpp_pypy.float_w_interleaved_warm", "ours_python.pypy_warm"),
    ("closing", "from", "PyPy", "JIT warm-up", "control", B7S, "ours_python.pypy_warm", "python"),
    # B8 e B9: Lua (LZ1, LZ2, L5), LuaJIT (L1-L5) e C# com os laços otimizados desde o início.
    ("closing", "from", "C++ -> Lua", "L1-L4", "ours_in_translation", B8I, "cpp_lua.lua_all_hoist", "lua"),
    ("closing", "from", "C++ -> Lua", "L1-L5", "ours_in_translation", B8I, "cpp_lua.lua_all_hoist_carry", "lua"),
    ("closing", "from", "C++ -> Lua", "L1-L4 (translation); index without -1 (ours)", "both", B8I, "cpp_lua.lua_all_hoist", "ours_lua.zero_based_index"),
    ("closing", "from", "C++ -> Lua", "L1-L5 (translation); index without -1 (ours)", "both", B8I, "cpp_lua.lua_all_hoist_carry", "ours_lua.zero_based_index"),
    ("closing", "from", "Lua", "Index without -1", "control", B8I, "ours_lua.zero_based_index", "lua"),
    ("closing", "from", "C++ -> LuaJIT", "Subsequence storage, B9", "ours_in_translation", B9C, "cpp_lua.luajit_subseq_array", "lua"),
    ("closing", "from", "C++ -> LuaJIT", "L1-L5", "ours_in_translation", B9C, "cpp_lua.luajit_all_hoist_carry", "lua"),
    ("closing", "from", "C++ -> LuaJIT", "Subsequence storage (translation); literal positions (ours), B9", "both", B9C, "cpp_lua.luajit_subseq_array", "ours_lua.literal_offsets"),
    ("closing", "from", "C++ -> LuaJIT", "L1-L5 (translation); literal positions (ours)", "both", B9C, "cpp_lua.luajit_all_hoist_carry", "ours_lua.literal_offsets"),
    ("closing", "from", "C++ -> C#", "Bounds checks + distance matrix, B9", "ours_in_translation", B9C, "cpp_csharp.unchecked_seq_flat_matrix", "csharp"),
    ("closing", "from", "C++ -> C#", "Bounds checks + distance matrix, loops optimized from the start in both", "both", B9C, "cpp_csharp.unchecked_seq_flat_matrix_quickjit_loops_off", "ours_csharp.quickjit_loops_off"),
    ("closing", "from", "C#", "Loops optimized from the start", "control", B9C, "ours_csharp.quickjit_loops_off", "csharp"),

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
    for round_, prog in (("base-20261120", "cplusplus"), (B4I, "lua")):
        with open(RUNS / round_ / "csv" / f"{prog}.csv", newline="") as f:
            for r in csv.DictReader(f):
                costs.setdefault(r["instance"].removesuffix(".tsp"), r["COST"])
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
        insts = instances(round_)
        ratios = [mean(a[i]) / mean(b[i]) for i in insts]
        cv = max(stdev(v) / mean(v) for t in (a, b) for v in (t[i] for i in insts))
        by_inst = dict(zip(insts, ratios))
        rows.append({
            "table": table, "direction": direction, "pair": pair, "factor": factor, "column": column,
            "round": round_, "numerator": num, "denominator": den,
            **{f"ratio_{i}": f"{by_inst[i]:.4f}" if i in by_inst else "" for i in ALL},
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

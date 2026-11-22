#!/usr/bin/env python3
"""
Benchmark de TODOS os projetos de traducao desta pasta.

Cada pasta `<origem>_<destino>` contem uma traducao do solver MLP para a
linguagem `<destino>`. Este script compara essa traducao com a implementacao
ORIGINAL DA MESMA LINGUAGEM que existe na pasta superior:

    claude/cpp_go     (Go traduzido do C++)  x  ../go     (Go original)
    claude/cpp_rust   (Rust traduzido)       x  ../rust   (Rust original)
    claude/java_cpp   (C++ traduzido)        x  ../cplusplus
    ...

Regra central: os metodos de compilacao/execucao sao IDENTICOS para os dois
lados. Cada toolchain aqui replica o `build.sh`/`makefile`/`run_*.sh` da
implementacao de referencia -- ver REFERENCIA_* abaixo, onde cada conjunto de
flags cita o arquivo de onde veio.

As referencias NAO sao compiladas dentro de `../<linguagem>/`: elas sao
copiadas para `_referencias/` e compiladas la, para nao mexer no repositorio
pai. Os solvers leem "../distance_matrix", entao a matriz e replicada tanto na
raiz desta pasta quanto na raiz de `_referencias/`.

Uso:
    python3 benchmark.py                      # tudo
    python3 benchmark.py -n 3                 # 3 repeticoes por medicao
    python3 benchmark.py -p cpp_rust cpp_go   # so alguns projetos
    python3 benchmark.py -i berlin52          # so algumas instancias
    python3 benchmark.py --no-build           # pula compilacao (e a copia das refs)
    python3 benchmark.py --incluir-original   # roda tambem o original da propria pasta

O Lua e medido nos dois runtimes da referencia (`lua` e `luajit`), gerando uma
linha por runtime -- ver `runtimes_extras` em PROJETOS.
"""

import argparse
import csv
import glob
import os
import queue
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

# =============================================================================
# CONFIGURACOES GERAIS
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if os.path.isdir(os.path.join(SCRIPT_DIR, "claude")):
    RAIZ = os.path.join(SCRIPT_DIR, "claude")
    PAI = SCRIPT_DIR
else:
    RAIZ = SCRIPT_DIR
    PAI = os.path.dirname(RAIZ)

INSTANCIAS = [
    "burma14", "ulysses16", "gr17", "gr21", "ulysses22", "gr24", "fri26",
    "bayg29", "bays29", "dantzig42", "swiss42", "att48", "gr48", "hk48",
    "eil51", "berlin52", "brazil58", "st70", "eil76", "pr76", "gr96",
    "rat99", "kroA100", "kroB100", "kroC100", "kroD100", "kroE100", "rd100",
    "eil101", "lin105", "pr107", "gr120", "pr124", "bier127", "ch130",
    "pr136", "gr137", "pr144", "ch150", "kroA150", "kroB150", "pr152",
    "u159", "si175", "brg180", "rat195", "d198", "kroA200", "kroB200",
    "gr202", "ts225", "tsp225", "pr226", "gr229", "gil262", "pr264",
    "a280", "pr299",
]
PASTA_INSTANCIAS = os.path.join(RAIZ, "mlp-instances")
PASTA_REFS = os.path.join(RAIZ, "_referencias")
PASTA_MATRIZES = os.path.join(RAIZ, "_matrizes")
PASTA_WORKERS = os.path.join(RAIZ, "_workers")
NIX_ENV = os.path.join(PAI, "nix-env.sh")

SAIDA_TXT = os.path.join(RAIZ, "resultados_mlp.txt")
SAIDA_CSV = os.path.join(RAIZ, "resultados_mlp.csv")

TIMEOUT_PADRAO = 7200  # segundos por execucao individual (igual a ../run_bm.py)

# artefatos de build / lixo que nao precisam ser copiados para _referencias
IGNORAR_NA_COPIA = shutil.ignore_patterns(
    "obj", "mod", "target", "bin", ".build", "__pycache__", "node_modules",
    "*.class", "*.o", "*.mod", "callgrind.out.*", ".git", "_workers", "_matrizes",
)

# =============================================================================
# FLAGS DE REFERENCIA (copiadas das implementacoes da pasta superior)
# =============================================================================

# ../cplusplus/makefile -> CCOPT
REF_CPP_FLAGS = [
    "-m64", "-Ofast", "-fPIC", "-fexceptions",
    "-DNDEBUG", "-DIL_STD", "-Winline", "-g",
]
# ../cplusplus/makefile usa -std=c++0x. Algumas traducoes usam construcoes
# posteriores; nesse caso o build e refeito com o proximo padrao da lista.
# As flags de otimizacao (que sao o que importa para a comparacao) nao mudam.
REF_CPP_STDS = ["-std=c++0x", "-std=c++17"]

# ../c/makefile -> COPT (MACRO = -DFLAT -DNDEBUG) e LOPT
REF_C_FLAGS = [
    "-std=c99", "-DFLAT", "-DNDEBUG",
    "-Wno-unused-result", "-Winline", "-g", "-Ofast",
]
REF_C_LINK = ["-lm"]

# ../fortran/makefile -> PREPROC + COPT
REF_F_FLAGS = ["-cpp", "-DDOUBLE", "-Winline", "-g", "-Ofast"]

# ../go/build.sh
REF_GO_FLAGS = ["-gcflags=-B"]

# ../java/build.sh e ../java/run_java.sh
REF_JAVAC_FLAGS = ["-O", "-Xlint:unchecked"]

# ../csharp/build.sh  (+ ../csharp/csharp.csproj nao liga ServerGC nem
# InvariantGlobalization, entao forcamos os mesmos defaults nos projetos que ligam)
REF_DOTNET_FLAGS = [
    "-c", "Release",
    "-p:ServerGarbageCollection=false",
    "-p:InvariantGlobalization=false",
]

# ../rust/build.sh + ../rust/Cargo.toml [profile.release]
REF_CARGO_FLAGS = [
    "--release",
    "--config", "profile.release.opt-level=3",
    "--config", "profile.release.codegen-units=1",
    "--config", "profile.release.debug=1",
    "--config", 'profile.release.panic="abort"',
]

# ../swift/build.sh
REF_SWIFT_BUILD_FLAGS = [
    "-c", "release", "-Xswiftc", "-Xfrontend", "-Xswiftc", "-Ounchecked",
]
# equivalente para os projetos sem Package.swift (arquivos soltos)
REF_SWIFTC_FLAGS = ["-O", "-Xfrontend", "-Ounchecked"]

# runtimes interpretados -- exatamente os run_*.sh da referencia
REF_NODE = ["node", "--no-concurrent-recompilation"]           # ../javascript/run_node.sh
REF_JULIA = ["julia", "--check-bounds=no", "--math-mode=fast"]  # ../julia/run_julia.sh
REF_LUA = ["lua"]                                               # ../lua/run_lua.sh
REF_LUAJIT = ["luajit"]                                         # ../lua/run_luajit.sh
REF_PYTHON = ["python3"]                                        # ../python/run_python3.sh
REF_PYPY = ["pypy3"]                                            # ../python/run_pypy.sh
REF_JAVA = ["java"]                                             # ../java/run_java.sh


def _cpp(bin_name, objdir="obj", srcs="*.cpp", extra=()):
    return ("cpp", {"srcs": srcs, "objdir": objdir, "bin": bin_name,
                    "extra": list(extra)})


# =============================================================================
# IMPLEMENTACOES DE REFERENCIA (pasta superior) -- a linha de base
# =============================================================================
# "src" e a pasta em ../ que sera copiada para _referencias/<src>.
# build/run replicam o build.sh e o run_*.sh daquela pasta.
REFERENCIAS = {
    "C": {
        "src": "c",
        "build": ("c", {"srcs": "*.c", "objdir": "obj", "bin": "solve"}),
        "run": ["./solve"],
    },
    "C++": {
        "src": "cplusplus",
        "build": _cpp("solve"),
        "run": ["./solve"],
    },
    "C#": {
        "src": "csharp",
        "build": ("dotnet", {"proj": "csharp.csproj"}),
        "run": ["./bin/Release/net8.0/csharp"],
    },
    "Fortran": {
        "src": "fortran",
        # ordem vinda das dependencias declaradas em ../fortran/makefile
        "build": ("fortran", {"srcs": ["types.f90", "assert.f90", "data.f90", "main.f90"],
                              "objdir": "obj", "moddir": "mod", "bin": "solve"}),
        "run": ["./solve"],
    },
    "Go": {
        "src": "go",
        "build": ("go", {"pkg": "main.go", "bin": "main"}),
        "run": ["./main"],
    },
    "Java": {
        "src": "java",
        "build": ("java", {"entry": "main.java"}),
        "run": REF_JAVA + ["Main"],
    },
    "JavaScript": {"src": "javascript", "build": None, "run": REF_NODE + ["main.js"]},
    "Julia":      {"src": "julia",      "build": None, "run": REF_JULIA + ["main.jl"]},
    "Lua":        {"src": "lua",        "build": None, "run": REF_LUA + ["main.lua"]},
    # mesmo fonte de "Lua", outro runtime -- igual ao par Python3/PyPy3 sobre ../python
    "LuaJIT":     {"src": "lua",        "build": None, "run": REF_LUAJIT + ["main.lua"]},
    "Python3":    {"src": "python",     "build": None, "run": REF_PYTHON + ["main.py"]},
    "PyPy3":      {"src": "python",     "build": None, "run": REF_PYPY + ["main_pypy.py"]},
    "Rust": {
        "src": "rust",
        "build": ("cargo", {"bin": "gils_rvnd", "features": ["flat"], "copy_to": "solve"}),
        "run": ["./solve"],
    },
    "Swift": {
        "src": "swift",
        "build": ("swiftpm", {}),
        "run": ["./.build/release/solve"],
    },
}


# =============================================================================
# CATALOGO DOS PROJETOS
# =============================================================================
# "traduzida" e o que sera comparado com REFERENCIAS[lang].
# "original" e a implementacao de origem que vive na mesma pasta; so roda com
# --incluir-original (ela nao faz parte da comparacao pedida).
#
# Todos os solvers leem "../distance_matrix", logo o cwd de execucao e SEMPRE
# a pasta do projeto.
PROJETOS = {
    # ---------------------------------------------------------------- C <-> C++
    "c_cpp": {
        "original": {"lang": "C",
                     "build": ("c", {"srcs": "*.c", "objdir": "obj", "bin": "solve"}),
                     "run": ["./solve"]},
        # fontes nao tem `#define FLAT` proprio (ao contrario dos cpp_*),
        # entao a macro vem pela linha de comando, como em ../c/makefile
        "traduzida": {"lang": "C++",
                      "build": _cpp("solve_cpp", objdir="obj_cpp", extra=["-DFLAT"]),
                      "run": ["./solve_cpp"]},
    },
    "cpp_c": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "C",
                      "build": ("c", {"srcs": "*.c", "objdir": "obj_c", "bin": "solve_c"}),
                      "run": ["./solve_c"]},
    },

    # ---------------------------------------------------- C++ -> outras linguagens
    "cpp_csharp": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        # detecta automaticamente se dotnet fica em cpp_csharp/csharp ou direto em cpp_csharp/
        "traduzida": {
            "lang": "C#",
            **({"subdir": "csharp", "build": ("dotnet", {"proj": "csharp.csproj"}), "run": ["./bin/Release/net8.0/csharp"]}
               if os.path.isdir(os.path.join(RAIZ, "cpp_csharp", "csharp"))
               else {"build": ("dotnet", {"proj": "mlp.csproj"}), "run": ["./bin/Release/net8.0/solve"]})
        },
    },
    "cpp_fortran": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "Fortran",
                      "build": ("fortran", {"srcs": ["mlp.f90"] if os.path.exists(os.path.join(RAIZ, "cpp_fortran", "mlp.f90"))
                                            else ["data_mod.f90", "read_data.f90", "main.f90"],
                                            "objdir": "obj_f", "moddir": "obj_f",
                                            "bin": "solve_f"}),
                      "run": ["./solve_f"]},
    },
    "cpp_go": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "Go",
                      "build": ("go", {"pkg": "main.go" if os.path.exists(os.path.join(RAIZ, "cpp_go", "main.go")) else "./go_src",
                                       "bin": "solve_go"}),
                      "run": ["./solve_go"]},
    },
    "cpp_java": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {
            "lang": "Java",
            **({"build": ("java", {"entry": "MLPSolver.java"}), "run": REF_JAVA + ["MLPSolver"]}
               if os.path.exists(os.path.join(RAIZ, "cpp_java", "MLPSolver.java"))
               else {"build": ("java", {"entry": "Main.java"}), "run": REF_JAVA + ["Main"]})
        },
    },
    "cpp_javascript": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "JavaScript", "build": None, "run": REF_NODE + ["main.js"]},
    },
    "cpp_julia": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "Julia", "build": None,
                      "run": REF_JULIA + (["MLP_Solver.jl"] if os.path.exists(os.path.join(RAIZ, "cpp_julia", "MLP_Solver.jl")) else ["main.jl"])},
    },
    "cpp_lua": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {
            "lang": "Lua", "build": None,
            "run": REF_LUA + (["solve.lua"] if os.path.exists(os.path.join(RAIZ, "cpp_lua", "solve.lua")) else ["main.lua"]),
            "runtimes_extras": [{"lang": "LuaJIT",
                                 "run": REF_LUAJIT + (["solve.lua"] if os.path.exists(os.path.join(RAIZ, "cpp_lua", "solve.lua")) else ["main.lua"])}]
        },
    },
    "cpp_pypy": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "PyPy3", "build": None,
                      "run": REF_PYPY + (["main.py"] if os.path.exists(os.path.join(RAIZ, "cpp_pypy", "main.py")) else ["mlp_solver.py"])},
    },
    "cpp_python": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {"lang": "Python3", "build": None, "run": REF_PYTHON + ["main.py"]},
    },
    "cpp_rust": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        "traduzida": {
            "lang": "Rust",
            **({"build": ("cargo", {"bin": "solve"}), "run": ["./target/release/solve"]}
               if os.path.exists(os.path.join(RAIZ, "cpp_rust", "target", "release", "solve"))
               or (os.path.exists(os.path.join(RAIZ, "cpp_rust", "Cargo.toml")) and 'name = "solve"' in open(os.path.join(RAIZ, "cpp_rust", "Cargo.toml"), errors='ignore').read())
               else {"build": ("cargo", {"bin": "cpp_to_rust"}), "run": ["./target/release/cpp_to_rust"]})
        },
    },
    "cpp_swift": {
        "original": {"lang": "C++", "build": _cpp("solve"), "run": ["./solve"]},
        # sem Package.swift: swiftc com o equivalente das flags de ../swift/build.sh
        "traduzida": {"lang": "Swift",
                      "build": ("swiftc", {"srcs": "*.swift", "bin": "solve_swift"}),
                      "run": ["./solve_swift"]},
    },

    # ---------------------------------------------------- outras linguagens -> C++
    "csharp_cpp": {
        "original": {"lang": "C#",
                     "build": ("dotnet", {"proj": "csharp.csproj"}),
                     "run": ["./bin/Release/net8.0/csharp"]},
        "traduzida": {"lang": "C++", "build": _cpp("mlp_cpp"), "run": ["./mlp_cpp"]},
    },
    "fortran_cpp": {
        "original": {"lang": "Fortran",
                     "build": ("fortran", {"srcs": ["types.f90", "assert.f90",
                                                    "data.f90", "main.f90"],
                                           "objdir": "obj", "moddir": "mod",
                                           "bin": "solve"}),
                     "run": ["./solve"]},
        # preproc.hpp escolhe float/double via DOUBLE, igual ../fortran/makefile
        "traduzida": {"lang": "C++",
                      "build": _cpp("solve_cpp", objdir="obj_cpp", extra=["-DDOUBLE"]),
                      "run": ["./solve_cpp"]},
    },
    "go_cpp": {
        "original": {"lang": "Go",
                     "build": ("go", {"pkg": "main.go", "bin": "main"}),
                     "run": ["./main"]},
        "traduzida": {"lang": "C++",
                      "build": _cpp("main_cpp", objdir="obj_cpp"),
                      "run": ["./main_cpp"]},
    },
    "java_cpp": {
        "original": {"lang": "Java",
                     "build": ("java", {"entry": "main.java"}),
                     "run": REF_JAVA + ["Main"]},
        "traduzida": {"lang": "C++",
                      "build": _cpp("main_cpp", objdir="obj_cpp"),
                      "run": ["./main_cpp"]},
    },
    "javascript_cpp": {
        "original": {"lang": "JavaScript", "build": None, "run": REF_NODE + ["main.js"]},
        "traduzida": {"lang": "C++", "build": _cpp("main"), "run": ["./main"]},
    },
    "julia_cpp": {
        "original": {"lang": "Julia", "build": None, "run": REF_JULIA + ["main.jl"]},
        "traduzida": {"lang": "C++", "build": _cpp("mlp"), "run": ["./mlp"]},
    },
    "lua_cpp": {
        "original": {"lang": "Lua", "build": None, "run": REF_LUA + ["main.lua"],
                     "runtimes_extras": [{"lang": "LuaJIT",
                                          "run": REF_LUAJIT + ["main.lua"]}]},
        "traduzida": {"lang": "C++", "build": _cpp("main"), "run": ["./main"]},
    },
    "python_cpp": {
        "original": {"lang": "Python3", "build": None, "run": REF_PYTHON + ["main.py"]},
        "traduzida": {"lang": "C++", "build": _cpp("main"), "run": ["./main"]},
    },
    "rust_cpp": {
        # ../rust/build.sh: cargo build --release --features flat && cp .../gils_rvnd solve
        "original": {"lang": "Rust",
                     "build": ("cargo", {"bin": "gils_rvnd", "features": ["flat"],
                                         "copy_to": "solve"}),
                     "run": ["./solve"]},
        "traduzida": {
            "lang": "C++",
            **({"build": _cpp("solve_cpp", objdir="obj_cpp", srcs="main.cpp"), "run": ["./solve_cpp"]}
               if os.path.exists(os.path.join(RAIZ, "rust_cpp", "main.cpp"))
               else {"build": _cpp("solve", objdir="cpp/obj", srcs="cpp/*.cpp", extra=["-DFLAT"]), "run": ["./cpp/solve"]})
        },
    },
    "swift_cpp": {
        "original": {"lang": "Swift", "build": ("swiftpm", {}),
                     "run": ["./.build/release/solve"]},
        "traduzida": {
            "lang": "C++",
            **({"build": _cpp("solve_cpp", objdir="obj_cpp", srcs="main.cpp"), "run": ["./solve_cpp"]}
               if os.path.exists(os.path.join(RAIZ, "swift_cpp", "main.cpp"))
               else {"build": _cpp("solve", objdir="src/obj", srcs="src/*.cpp"), "run": ["./solve"]})
        },
    },
}


# =============================================================================
# AMBIENTE (nix)
# =============================================================================
def carregar_ambiente():
    """Resolve o ambiente do ../nix-env.sh UMA vez e reusa em todos os comandos.

    Assim nenhum `source` entra no caminho medido de cada execucao.
    """
    base = dict(os.environ)
    if not os.path.isfile(NIX_ENV):
        print(f"[aviso] {NIX_ENV} nao encontrado; usando o PATH atual.")
        env = base
    else:
        res = subprocess.run(
            ["bash", "-c", f'source "{NIX_ENV}" >/dev/null 2>&1 && env -0'],
            capture_output=True,
        )
        if res.returncode != 0 or not res.stdout:
            print("[aviso] falha ao carregar nix-env.sh; usando o PATH atual.")
            env = base
        else:
            env = {}
            for entry in res.stdout.decode("utf-8", "replace").split("\0"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    if k not in ("PWD", "OLDPWD", "SHLVL", "_"):
                        env[k] = v
            env = env or base

    # Locale unico e deterministico para todas as linguagens: tira do caminho a
    # formatacao numerica dependente de locale e evita que runtimes sem o
    # locale do sistema (pypy3 x pt_BR.UTF-8, aqui) abortem na inicializacao.
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return env


ENV = {}  # preenchido no main()


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, env=ENV, capture_output=True, text=True)


def _falha(res, cmd):
    return "$ {}\n{}\n{}".format(" ".join(cmd), res.stdout.strip(), res.stderr.strip())


# =============================================================================
# COMPILADORES (um por toolchain, espelhando a referencia)
# =============================================================================
def build_cpp(cwd, spec):
    """Igual a ../cplusplus/makefile: um .o por .cpp e depois o link."""
    srcs = spec["srcs"]
    srcs = sorted(glob.glob(os.path.join(cwd, srcs))) if isinstance(srcs, str) \
        else [os.path.join(cwd, s) for s in srcs]
    if not srcs:
        return False, "nenhum fonte .cpp encontrado", None

    objdir = os.path.join(cwd, spec["objdir"])
    extra = spec.get("extra", [])
    ultimo = ""

    for std in REF_CPP_STDS:
        shutil.rmtree(objdir, ignore_errors=True)
        os.makedirs(objdir, exist_ok=True)
        flags = REF_CPP_FLAGS + [std] + extra
        objs, erro = [], None

        for src in srcs:
            obj = os.path.join(objdir, os.path.basename(src)[:-4] + ".o")
            cmd = ["g++"] + flags + ["-c", src, "-o", obj]
            res = _run(cmd, cwd)
            if res.returncode != 0:
                erro = _falha(res, cmd)
                break
            objs.append(obj)

        if erro is None:
            binario = os.path.join(cwd, spec["bin"])
            cmd = ["g++", "-m64"] + objs + ["-o", binario] + flags
            res = _run(cmd, cwd)
            if res.returncode == 0:
                nota = None if std == REF_CPP_STDS[0] else f"compilado com {std}"
                return True, None, nota
            erro = _falha(res, cmd)

        ultimo = erro

    return False, ultimo, None


def build_c(cwd, spec):
    """Igual a ../c/makefile: um .o por .c, link com -lm."""
    srcs = sorted(glob.glob(os.path.join(cwd, spec["srcs"])))
    if not srcs:
        return False, "nenhum fonte .c encontrado", None

    objdir = os.path.join(cwd, spec["objdir"])
    shutil.rmtree(objdir, ignore_errors=True)
    os.makedirs(objdir, exist_ok=True)

    objs = []
    for src in srcs:
        obj = os.path.join(objdir, os.path.basename(src)[:-2] + ".o")
        cmd = ["gcc"] + REF_C_FLAGS + ["-c", src, "-o", obj]
        res = _run(cmd, cwd)
        if res.returncode != 0:
            return False, _falha(res, cmd), None
        objs.append(obj)

    binario = os.path.join(cwd, spec["bin"])
    cmd = ["gcc"] + objs + REF_C_LINK + ["-o", binario]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_fortran(cwd, spec):
    """Igual a ../fortran/makefile: -cpp -DDOUBLE -Ofast, modulos em moddir."""
    objdir = os.path.join(cwd, spec["objdir"])
    moddir = os.path.join(cwd, spec["moddir"])
    shutil.rmtree(objdir, ignore_errors=True)
    os.makedirs(objdir, exist_ok=True)
    os.makedirs(moddir, exist_ok=True)

    objs = []
    for nome in spec["srcs"]:  # ordem importa: modulos antes de quem os usa
        src = os.path.join(cwd, nome)
        obj = os.path.join(objdir, nome[:-4] + ".o")
        cmd = (["gfortran"] + REF_F_FLAGS
               + ["-c", src, "-o", obj, "-I" + moddir, "-J" + moddir])
        # compila DENTRO de moddir: o gfortran procura .mod tambem no diretorio
        # corrente, e um .mod velho na pasta do projeto (gerado por outra versao
        # do compilador) tem precedencia e aborta a compilacao
        res = _run(cmd, moddir)
        if res.returncode != 0:
            return False, _falha(res, cmd), None
        objs.append(obj)

    binario = os.path.join(cwd, spec["bin"])
    cmd = ["gfortran"] + objs + ["-o", binario, "-g"]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_go(cwd, spec):
    """Igual a ../go/build.sh: go build -gcflags=-B."""
    cmd = ["go", "build"] + REF_GO_FLAGS + ["-o", spec["bin"], spec["pkg"]]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_java(cwd, spec):
    """Igual a ../java/build.sh: rm *.class && javac -O <entry> -Xlint:unchecked."""
    for c in glob.glob(os.path.join(cwd, "*.class")):
        os.remove(c)
    cmd = ["javac"] + REF_JAVAC_FLAGS + [spec["entry"]]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_dotnet(cwd, spec):
    """Igual a ../csharp/build.sh: dotnet build -c Release."""
    cmd = ["dotnet", "build"] + REF_DOTNET_FLAGS + [spec["proj"]]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_cargo(cwd, spec):
    """Igual a ../rust/build.sh: rm -rf target && cargo build --release [--features flat]."""
    shutil.rmtree(os.path.join(cwd, "target"), ignore_errors=True)
    cmd = ["cargo", "build"] + REF_CARGO_FLAGS
    if spec.get("features"):
        cmd += ["--features", ",".join(spec["features"])]
    res = _run(cmd, cwd)
    if res.returncode != 0:
        return False, _falha(res, cmd), None
    if spec.get("copy_to"):
        shutil.copy2(os.path.join(cwd, "target", "release", spec["bin"]),
                     os.path.join(cwd, spec["copy_to"]))
    return True, None, None


def build_swiftpm(cwd, spec):
    """Igual a ../swift/build.sh."""
    cmd = ["swift", "build"] + REF_SWIFT_BUILD_FLAGS
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


def build_swiftc(cwd, spec):
    """Projeto sem Package.swift: swiftc com o equivalente de ../swift/build.sh."""
    srcs = sorted(glob.glob(os.path.join(cwd, spec["srcs"])))
    if not srcs:
        return False, "nenhum fonte .swift encontrado", None
    cmd = ["swiftc"] + REF_SWIFTC_FLAGS + srcs + ["-o", spec["bin"]]
    res = _run(cmd, cwd)
    return (True, None, None) if res.returncode == 0 else (False, _falha(res, cmd), None)


TOOLCHAINS = {
    "cpp": build_cpp,
    "c": build_c,
    "fortran": build_fortran,
    "go": build_go,
    "java": build_java,
    "dotnet": build_dotnet,
    "cargo": build_cargo,
    "swiftpm": build_swiftpm,
    "swiftc": build_swiftc,
}


# =============================================================================
# ALVOS
# =============================================================================
# Um "alvo" e qualquer coisa que precise ser compilada e/ou executada.
#   chave  : identificador unico
#   rotulo : como aparece no log
#   cwd    : diretorio de trabalho (a matriz e sempre lida em ../distance_matrix)
def alvo_referencia(lang):
    ref = REFERENCIAS[lang]
    return {
        "chave": f"ref:{lang}",
        "rotulo": f"../{ref['src']} [{lang}]",
        "cwd": os.path.join(PASTA_REFS, ref["src"]),
        "build": ref["build"], "run": ref["run"],
        "papel": "referencia", "lang": lang, "projeto": f"../{ref['src']}",
    }


# os dois runtimes que executam o mesmo fonte Lua
RUNTIMES_LUA = ("Lua", "LuaJIT")
ESCOLHAS_LUA = {"ambos": RUNTIMES_LUA, "lua": ("Lua",), "luajit": ("LuaJIT",)}


def alvo_projeto(proj, papel, lua_runtime="ambos"):
    """Devolve uma LISTA: uma variante pode ser medida em mais de um runtime.

    E o caso do Lua, que roda no interpretador `lua` e no JIT `luajit` a partir
    do mesmo fonte. `lua_runtime` filtra qual(is) medir. O build (quando existe)
    e feito uma unica vez, no primeiro alvo que sobreviver ao filtro.
    """
    var = PROJETOS[proj][papel]
    candidatos = [{"lang": var["lang"], "run": var["run"], "chave": f"{proj}:{papel}"}]
    for extra in var.get("runtimes_extras", []):
        candidatos.append({"lang": extra["lang"], "run": extra["run"],
                           "chave": f"{proj}:{papel}:{extra['lang']}"})

    permitidos = ESCOLHAS_LUA[lua_runtime]
    candidatos = [c for c in candidatos
                  if c["lang"] not in RUNTIMES_LUA or c["lang"] in permitidos]

    alvos = []
    for i, c in enumerate(candidatos):
        alvos.append({
            "chave": c["chave"],
            "rotulo": f"{proj} [{c['lang']} / {papel}]",
            # "subdir": projeto cujo build/execucao acontece numa subpasta
            "cwd": os.path.join(RAIZ, proj, var.get("subdir", "")),
            # so o primeiro alvo compila; os demais reusam o mesmo artefato
            "build": var["build"] if i == 0 else None,
            "run": c["run"],
            "papel": papel, "lang": c["lang"], "projeto": proj,
        })
    return alvos


def montar_alvos(projetos, incluir_original, lua_runtime="ambos"):
    """Traducoes + as referencias que elas exigem (cada referencia so uma vez)."""
    alvos, vistos = [], set()

    for proj in projetos:
        for a in alvo_projeto(proj, "traduzida", lua_runtime):
            alvos.append(a)
            vistos.add(a["chave"])

            lang = a["lang"]
            if lang not in REFERENCIAS:
                sys.exit(f"[erro] sem referencia em ../ para a linguagem '{lang}' ({proj})")
            r = alvo_referencia(lang)
            if r["chave"] not in vistos:
                alvos.append(r)
                vistos.add(r["chave"])

        if incluir_original:
            for o in alvo_projeto(proj, "original", lua_runtime):
                alvos.append(o)
                vistos.add(o["chave"])

    return alvos


# =============================================================================
# FASE 0 -- COPIA DAS REFERENCIAS
# =============================================================================
def sincronizar_referencias(alvos):
    """Copia ../<lang>/ para _referencias/<lang>/ (sem artefatos de build).

    Compilar a copia, e nao o original, deixa o repositorio pai intacto.
    """
    dirs = sorted({REFERENCIAS[a["lang"]]["src"] for a in alvos
                   if a["papel"] == "referencia"})
    if not dirs:
        return

    print(f"Fase 0: copiando {len(dirs)} referencia(s) de ../ para _referencias/\n")
    os.makedirs(PASTA_REFS, exist_ok=True)

    for d in dirs:
        origem = os.path.join(PAI, d)
        destino = os.path.join(PASTA_REFS, d)
        if not os.path.isdir(origem):
            print(f"  [FALHA] ../{d} nao existe")
            continue
        shutil.rmtree(destino, ignore_errors=True)
        shutil.copytree(origem, destino, ignore=IGNORAR_NA_COPIA, symlinks=True)
        print(f"  [ok]   ../{d}  ->  _referencias/{d}")
    print()


# =============================================================================
# FASE 1 -- COMPILACAO
# =============================================================================
def compilar_tudo(alvos, jobs):
    print("Fase 1: compilacao (mesmas flags dos dois lados)\n")
    falhas = {}

    def tarefa(alvo):
        if alvo["build"] is None:
            return alvo, True, None, "interpretada"
        toolchain, spec = alvo["build"]
        ok, err, nota = TOOLCHAINS[toolchain](alvo["cwd"], spec)
        return alvo, ok, err, nota

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for alvo, ok, err, nota in ex.map(tarefa, alvos):
            if ok:
                print(f"  [ok]   {alvo['rotulo']}" + (f"  ({nota})" if nota else ""))
                alvo["nota_build"] = nota
            else:
                print(f"  [FALHA] {alvo['rotulo']}")
                falhas[alvo["rotulo"]] = err
                alvo["quebrada"] = True

    if falhas:
        print("\n  --- detalhes das falhas de compilacao ---")
        for rotulo, err in falhas.items():
            print(f"\n  ### {rotulo}\n{err}\n")
    return falhas


# =============================================================================
# FASE 2 -- MATRIZES DE DISTANCIA E WORKERS
# =============================================================================
def preparar_matrizes(instancias):
    """Pre-carrega as matrizes de distancia para _matrizes/<inst>.matrix."""
    os.makedirs(PASTA_MATRIZES, exist_ok=True)
    faltando = [i for i in instancias if not os.path.exists(os.path.join(PASTA_MATRIZES, f"{i}.matrix"))]
    if faltando:
        print(f"Fase 2: pre-gerando {len(faltando)} matriz(es) de distancia em {PASTA_MATRIZES}...\n")
        for inst in faltando:
            cmd = ["./load", f"{inst}.tsp"]
            res = subprocess.run(cmd, cwd=PASTA_INSTANCIAS, env=ENV, capture_output=True, text=True)
            matriz = os.path.join(RAIZ, "distance_matrix")
            if res.returncode != 0 or not os.path.exists(matriz):
                print(f"  [erro] load falhou para {inst}: {res.stdout.strip()} {res.stderr.strip()}")
                continue
            dest = os.path.join(PASTA_MATRIZES, f"{inst}.matrix")
            shutil.copy2(matriz, dest)
            if os.path.exists(matriz):
                os.remove(matriz)


def preparar_worker_dir(worker_dir, alvos):
    """Cria a estrutura de diretorios do worker com symlinks para os artefatos de cada alvo."""
    os.makedirs(worker_dir, exist_ok=True)
    refs_dir = os.path.join(worker_dir, "_referencias")
    os.makedirs(refs_dir, exist_ok=True)
    link_ref = os.path.join(refs_dir, "distance_matrix")
    if not os.path.lexists(link_ref):
        os.symlink("../distance_matrix", link_ref)

    for alvo in alvos:
        rel_cwd = os.path.relpath(alvo["cwd"], RAIZ)
        w_cwd = os.path.join(worker_dir, rel_cwd)
        os.makedirs(w_cwd, exist_ok=True)
        if os.path.exists(alvo["cwd"]):
            for item in os.listdir(alvo["cwd"]):
                src = os.path.join(alvo["cwd"], item)
                dst = os.path.join(w_cwd, item)
                if not os.path.lexists(dst):
                    os.symlink(src, dst)


def carregar_instancia_para_worker(worker_dir, inst):
    """Atualiza worker_dir/distance_matrix para apontar para a matriz de inst."""
    matriz_src = os.path.join(PASTA_MATRIZES, f"{inst}.matrix")
    matriz_link = os.path.join(worker_dir, "distance_matrix")
    if os.path.islink(matriz_link):
        try:
            if os.readlink(matriz_link) == matriz_src:
                return True
        except OSError:
            pass
        os.remove(matriz_link)
    elif os.path.exists(matriz_link):
        os.remove(matriz_link)

    if not os.path.exists(matriz_src):
        return False
    os.symlink(matriz_src, matriz_link)
    return True


# =============================================================================
# FASE 3 -- EXECUCAO
# =============================================================================
RE_COST = re.compile(r"COST\s*:\s*([-\d.,eE+]+)")
RE_TIME = re.compile(r"TIME\s*:\s*([-\d.,eE+]+)")


def _extrair(padrao, texto):
    m = padrao.search(texto)
    if not m:
        return None
    try:
        return float(m.group(1).replace(",", "."))
    except ValueError:
        return None


def executar(alvo, timeout, cwd=None):
    if cwd is None:
        cwd = alvo["cwd"]
    inicio = time.perf_counter()
    try:
        res = subprocess.run(alvo["run"], cwd=cwd, env=ENV, capture_output=True,
                             text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": f"timeout ({timeout}s)", "wall": None, "time": None, "cost": None}
    except Exception as e:
        return {"status": f"erro: {e}", "wall": None, "time": None, "cost": None}

    wall = time.perf_counter() - inicio
    saida = res.stdout + res.stderr
    if res.returncode != 0:
        return {"status": f"erro ({res.returncode})", "wall": wall,
                "time": None, "cost": None, "saida": saida}
    return {"status": "ok", "wall": wall,
            "time": _extrair(RE_TIME, saida), "cost": _extrair(RE_COST, saida)}


def rodar_tudo(alvos, instancias, repeticoes, timeout, jobs=1):
    preparar_matrizes(instancias)

    print(f"\nFase 3: execucao paralela com {jobs} worker(s)\n")

    # Prepara os diretorios de cada worker
    os.makedirs(PASTA_WORKERS, exist_ok=True)
    worker_dirs = []
    worker_queue = queue.Queue()
    for w_id in range(jobs):
        w_dir = os.path.join(PASTA_WORKERS, f"worker_{w_id}")
        preparar_worker_dir(w_dir, alvos)
        worker_dirs.append(w_dir)
        worker_queue.put(w_id)

    # Lista de tarefas: (inst, alvo)
    tarefas = []
    for inst in instancias:
        mat_path = os.path.join(PASTA_MATRIZES, f"{inst}.matrix")
        if not os.path.exists(mat_path):
            continue
        for alvo in alvos:
            if not alvo.get("quebrada"):
                tarefas.append((inst, alvo))

    total_tarefas = len(tarefas)
    concluidas = 0
    lock = threading.Lock()
    registros = []

    campos = ["instancia", "projeto", "linguagem", "papel", "custo",
              "tempo_cpu_medio", "tempo_cpu_min", "tempo_wall_medio",
              "tempo_wall_min", "repeticoes", "status", "nota"]

    # Inicializa CSV com cabecalho (escrita incremental thread-safe)
    with lock:
        with open(SAIDA_CSV, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=campos)
            w.writeheader()

    def processar_tarefa(item):
        nonlocal concluidas
        inst, alvo = item
        w_id = worker_queue.get()
        try:
            w_dir = worker_dirs[w_id]
            if not carregar_instancia_para_worker(w_dir, inst):
                return None

            rel_cwd = os.path.relpath(alvo["cwd"], RAIZ)
            w_cwd = os.path.join(w_dir, rel_cwd)

            medicoes, status, custo = [], "ok", None
            for _ in range(repeticoes):
                r = executar(alvo, timeout, cwd=w_cwd)
                if r["status"] != "ok":
                    status = r["status"]
                    break
                custo = r["cost"] if custo is None else custo
                medicoes.append(r)

            reg = {
                "instancia": inst, "projeto": alvo["projeto"], "linguagem": alvo["lang"],
                "papel": alvo["papel"], "custo": custo,
                "tempo_cpu_medio": None, "tempo_cpu_min": None,
                "tempo_wall_medio": None, "tempo_wall_min": None,
                "repeticoes": 0, "status": status,
                "nota": alvo.get("nota_build") or "",
            }
            if status == "ok" and medicoes:
                walls = [m["wall"] for m in medicoes]
                times = [m["time"] for m in medicoes if m["time"] is not None]
                reg.update({
                    "tempo_cpu_medio": statistics.mean(times) if times else None,
                    "tempo_cpu_min": min(times) if times else None,
                    "tempo_wall_medio": statistics.mean(walls),
                    "tempo_wall_min": min(walls),
                    "repeticoes": len(medicoes),
                })

            with lock:
                concluidas += 1
                registros.append(reg)
                with open(SAIDA_CSV, "a", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=campos)
                    w.writerow(reg)
                if SCRIPT_DIR != RAIZ:
                    try:
                        shutil.copy2(SAIDA_CSV, os.path.join(SCRIPT_DIR, "resultados_mlp.csv"))
                    except Exception:
                        pass

                t = reg["tempo_cpu_medio"] if reg["tempo_cpu_medio"] is not None \
                    else reg["tempo_wall_medio"]
                tempo_str = f"{t:9.4f}s" if t is not None else "     -   "
                custo_str = f"COST {custo}" if custo is not None else ""
                progresso = f"[{concluidas}/{total_tarefas}]"
                status_str = f"({status})" if status != "ok" else ""
                print(f"  {progresso:<10} [w{w_id:<2}] {inst:<10} {alvo['rotulo'].ljust(44)} {tempo_str}  {custo_str} {status_str}")

            return reg
        finally:
            worker_queue.put(w_id)

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        list(ex.map(processar_tarefa, tarefas))

    try:
        shutil.rmtree(PASTA_WORKERS)
    except Exception:
        pass

    return registros


# =============================================================================
# RELATORIO
# =============================================================================
def _fmt(v, casas=4):
    return "-" if v is None else f"{v:.{casas}f}"


def _tempo(reg):
    if reg is None:
        return None
    return reg["tempo_cpu_medio"] if reg["tempo_cpu_medio"] is not None \
        else reg["tempo_wall_medio"]


def gerar_relatorio(registros, instancias, projetos, incluir_original):
    with open(SAIDA_CSV, "w", newline="", encoding="utf-8") as f:
        campos = ["instancia", "projeto", "linguagem", "papel", "custo",
                  "tempo_cpu_medio", "tempo_cpu_min", "tempo_wall_medio",
                  "tempo_wall_min", "repeticoes", "status", "nota"]
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(registros)

    def achar(inst, projeto, papel, lang=None):
        # a linguagem faz parte da chave: ../python e referencia de Python3 E de
        # PyPy3, entao so a pasta nao identifica o registro
        return next((r for r in registros if r["instancia"] == inst
                     and r["projeto"] == projeto and r["papel"] == papel
                     and (lang is None or r["linguagem"] == lang)), None)

    def achar_ref(inst, lang):
        return achar(inst, f"../{REFERENCIAS[lang]['src']}", "referencia", lang)

    def achar_todas(inst, projeto, papel):
        # um projeto pode ter mais de um registro por papel quando a mesma
        # variante e medida em varios runtimes (Lua x LuaJIT)
        return [r for r in registros if r["instancia"] == inst
                and r["projeto"] == projeto and r["papel"] == papel]

    largura = 96 if incluir_original else 82
    L = ["=" * largura,
         "RESULTADOS MLP -- traducao x implementacao original da MESMA linguagem (pasta ../)",
         "=" * largura, ""]

    for inst in instancias:
        if not any(r["instancia"] == inst for r in registros):
            continue
        L.append(f"Instancia: {inst}")
        L.append("-" * largura)
        cab = (f"{'Projeto':<18}{'Linguagem':<13}{'Referencia':<16}"
               f"{'Traduzida':>12}{'Original ../':>14}{'Razao':>9}")
        if incluir_original:
            cab += f"{'Origem/pasta':>14}"
        L.append(cab)
        L.append("-" * largura)

        for proj in projetos:
            # uma linha por runtime medido (Lua e LuaJIT geram duas)
            for trad in achar_todas(inst, proj, "traduzida"):
                ref = achar_ref(inst, trad["linguagem"])
                tt, tr = _tempo(trad), _tempo(ref)
                razao = f"{tt / tr:.2f}x" if (tt and tr) else "-"
                col_t = _fmt(tt) if tt else trad["status"][:11]
                col_r = _fmt(tr) if tr else (ref["status"][:13] if ref else "-")
                linha = (f"{proj:<18}{trad['linguagem']:<13}"
                         f"{'../' + REFERENCIAS[trad['linguagem']]['src']:<16}"
                         f"{col_t:>12}{col_r:>14}{razao:>9}")
                if incluir_original:
                    # casa pela linguagem quando existe; senao usa o primeiro
                    originais = achar_todas(inst, proj, "original")
                    orig = next((o for o in originais
                                 if o["linguagem"] == trad["linguagem"]),
                                originais[0] if originais else None)
                    linha += f"{_fmt(_tempo(orig)):>14}"
                L.append(linha)
        L.append("")
        L.append("  Razao = tempo da traducao / tempo do original da mesma linguagem.")
        L.append("  < 1.00x: a traducao ficou mais rapida que o original.")
        L.append("")

    L += ["=" * largura,
          "CUSTOS (traducao x original da mesma linguagem -- devem bater, mesma semente .rnd)",
          "=" * largura]
    for inst in instancias:
        if not any(r["instancia"] == inst for r in registros):
            continue
        L.append(f"\nInstancia: {inst}")
        for proj in projetos:
            for trad in achar_todas(inst, proj, "traduzida"):
                ref = achar_ref(inst, trad["linguagem"])
                ct = trad["custo"]
                cr = ref["custo"] if ref else None
                marca = "OK " if (ct is not None and cr is not None
                                  and abs(ct - cr) < 1e-6) else "!! "
                rot = f"{proj} ({trad['linguagem']})" if len(
                    achar_todas(inst, proj, "traduzida")) > 1 else proj
                L.append(f"  {marca}{rot:<26} traduzida={ct}  "
                         f"original=({trad['linguagem']}) {cr}")

    texto = "\n".join(L) + "\n"
    with open(SAIDA_TXT, "w", encoding="utf-8") as f:
        f.write(texto)
    if SCRIPT_DIR != RAIZ:
        try:
            shutil.copy2(SAIDA_TXT, os.path.join(SCRIPT_DIR, "resultados_mlp.txt"))
            shutil.copy2(SAIDA_CSV, os.path.join(SCRIPT_DIR, "resultados_mlp.csv"))
        except Exception:
            pass
    print(texto)
    print(f"[ok] relatorio salvo em {SAIDA_TXT} e {SAIDA_CSV}")


# =============================================================================
def main():
    global ENV

    ap = argparse.ArgumentParser(
        description="Compara cada traducao com a implementacao original da mesma linguagem")
    ap.add_argument("-n", "--runs", type=int, default=1,
                    help="repeticoes por (alvo, instancia)")
    ap.add_argument("-p", "--projects", nargs="+", default=None,
                    help="subconjunto de projetos (padrao: todos)")
    ap.add_argument("-i", "--instances", nargs="+", default=None,
                    help=f"subconjunto de instancias (padrao: {' '.join(INSTANCIAS)})")
    ap.add_argument("--incluir-original", action="store_true",
                    help="roda tambem a implementacao de origem que vive na propria pasta")
    ap.add_argument("--lua-runtime", choices=["ambos", "lua", "luajit"], default="ambos",
                    help="qual runtime Lua medir (padrao: ambos). 'lua' usa so o "
                         "interpretador, 'luajit' so o JIT -- vale para a traducao "
                         "e para a referencia ../lua")
    ap.add_argument("--no-build", action="store_true",
                    help="pula a copia das referencias e a compilacao")
    ap.add_argument("--build-only", action="store_true", help="so compila, nao executa")
    ap.add_argument("--timeout", type=int, default=TIMEOUT_PADRAO,
                    help="timeout por execucao, em segundos")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2),
                    help="grau de paralelismo para compilacao e execucao dos benchmarks (workers paralelos)")
    args = ap.parse_args()

    projetos = args.projects or sorted(PROJETOS)
    desconhecidos = [p for p in projetos if p not in PROJETOS]
    if desconhecidos:
        sys.exit(f"[erro] projeto(s) desconhecido(s): {', '.join(desconhecidos)}")

    instancias = args.instances or INSTANCIAS
    faltando = [i for i in instancias
                if not os.path.isfile(os.path.join(PASTA_INSTANCIAS, f"{i}.tsp"))]
    if faltando:
        sys.exit(f"[erro] instancia(s) nao encontrada(s): {', '.join(faltando)}")

    ENV = carregar_ambiente()
    alvos = montar_alvos(projetos, args.incluir_original, args.lua_runtime)

    if not args.no_build:
        sincronizar_referencias(alvos)
        compilar_tudo(alvos, args.jobs)
    if args.build_only:
        return

    registros = rodar_tudo(alvos, instancias, args.runs, args.timeout, jobs=args.jobs)

    # Ordena registros pela ordem das instancias e dos projetos para o relatorio final
    def chave_ordem(r):
        i_idx = instancias.index(r["instancia"]) if r["instancia"] in instancias else 9999
        p_idx = projetos.index(r["projeto"]) if r["projeto"] in projetos else 9999
        return (i_idx, p_idx, r["linguagem"], r["papel"])

    registros_ordenados = sorted(registros, key=chave_ordem)
    gerar_relatorio(registros_ordenados, instancias, projetos, args.incluir_original)


if __name__ == "__main__":
    main()

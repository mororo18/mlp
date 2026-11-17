# minimum latency problem

## setup
The steps to set up the enviroment are described [here](https://github.com/mororo18/mlp/blob/main/SETUP.md).

## how to build (if necessary)

example:
```
cd java
./build.sh
```

## run

example:
```
cd mlp-instances/
./load att48.tsp
cd ../java
./run_java.sh
```

## Build & Run Scripts Status

| Language | build.sh | run_*.sh | `-v`/`--verbose` |
|----------|----------|----------|------------------|
| c | ✓ | ✓ run_c.sh | ✓ |
| c_asm | ✓ | ✓ run_c.sh | — (not planned) |
| cplusplus | ✓ | ✓ run_cpp.sh | ✓ |
| cppOOP | ✓ | ✓ run_cpp-OOP.sh | ✓ |
| csharp | ✓ | ✓ run_dotnet.sh, run_mcs.sh | ✓ (via `run_dotnet.sh`; `run_mcs.sh`'s mono build is currently broken, see `TODO.md`) |
| fortran | ✓ | ✓ run_fortran.sh | ✓ |
| go | ✓ | ✓ run_golang.sh | ✓ |
| java | ✓ | ✓ run_java.sh | ✓ |
| javascript | ✗ | ✓ run_node.sh | ✓ |
| julia | ✗ | ✓ run_julia.sh | ✓ |
| lua | ✗ | ✓ run_lua.sh, run_luajit.sh | ✗ blocked, run scripts point to a broken entry point (see `TODO.md`) |
| octave | ✗ | ✓ run_octave.sh, run_matlab.sh | ✗ blocked, implementation never runs to completion (see `TODO.md`) |
| python | ✗ | ✓ run_python3.sh, run_pypy.sh | ✓ |
| rust | ✓ | ✓ run_rust.sh | ✓ |

**Note:** Interpreted languages (javascript, julia, lua, octave, python) don't require build.sh.

## Progress output (`-v`/`--verbose`)

Every solver accepts an optional `-v`/`--verbose` flag (default off) that
prints per-iteration GILS-RVND progress (construction cost, best-so-far
cost, current solution). Without it, only the final `COST`/`TIME` (and
similar per-language final metrics) are printed. See the table above for
per-language status.

```
cd c
./run_c.sh --verbose
```

## Tool Versions

| Tool | Version |
|------|---------|
| gcc | 14.2.0 |
| g++ | 14.2.0 |
| gfortran | 14.2.0 |
| go | 1.22.2 |
| rust | 1.66.1 (via rust-toolchain) |
| java (OpenJDK) | 21 |
| node | 20.16.0 |
| julia | 1.12 (via juliaup) |
| python3 | 3.12.7 |
| pypy3 | 7.3.15 |
| lua5.3 | 5.3.6 |
| luajit | 2.1.1703358377 |
| octave | 8.4.0 |
| dotnet | 8.0.125 |
| valgrind | 3.23.0 |

**Note:** Run `./setup.sh` to install all required tools.

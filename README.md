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

| Language | build.sh | run_*.sh |
|----------|----------|----------|
| c | ✓ | ✓ run_c.sh |
| c_asm | ✓ | ✓ run_c.sh |
| cplusplus | ✓ | ✓ run_cpp.sh |
| cppOOP | ✓ | ✓ run_cpp-OOP.sh |
| csharp | ✓ | ✓ run_dotnet.sh, run_mcs.sh |
| fortran | ✓ | ✓ run_fortran.sh |
| go | ✓ | ✓ run_golang.sh |
| java | ✓ | ✓ run_java.sh |
| javascript | ✗ | ✓ run_node.sh |
| julia | ✗ | ✓ run_julia.sh |
| lua | ✗ | ✓ run_lua.sh, run_luajit.sh |
| octave | ✗ | ✓ run_octave.sh, run_matlab.sh |
| python | ✗ | ✓ run_python3.sh, run_pypy.sh |
| rust | ✓ | ✓ run_rust.sh |

**Note:** Interpreted languages (javascript, julia, lua, octave, python) don't require build.sh.

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

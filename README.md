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

# MLG - Minimum Latency Problem

## Project Overview
- **MLP** = Minimum Latency Problem
- Algorithm: **GILS-RVND** (Greedy Inter-roulette Selection + Random Variable Neighborhood Descent)
- Goal: Benchmark the same algorithm across 15+ programming languages

## Repository Structure
```
mlp/
├── c/              - C implementation
├── c_asm/          - C + Assembly
├── cplusplus/      - C++ implementation
├── cppOOP/         - C++ OOP style
├── csharp/         - C# implementation
├── fortran/        - Fortran implementation
├── go/             - Go implementation
├── java/           - Java implementation
├── javascript/     - JavaScript/Node.js
├── julia/          - Julia implementation
├── lua/            - Lua implementation
├── octave/         - Octave/MATLAB
├── python/         - Python/PyPy implementations
├── rust/           - Rust implementation
├── mlp-instances/  - Test data (git submodule)
│   └── loader/     - C++ utility to load TSP instances
├── run_bm.py       - Benchmark runner
└── manager_bm.py   - Benchmark manager
```

## Setup
1. Initialize submodule:
   ```bash
   git submodule update --init
   ```
2. Install required tools (see **SETUP.md**)

## Build & Run Conventions
- **Build**: Each language directory has a `build.sh` script
- **Run**: Each language directory has a `run_*.sh` script
- **Pattern**: Load instance → Run solver
- Note: Interpreted languages (lua, python, octave, julia) may skip build.sh

Example:
```bash
cd mlp-instances/
./loader/load att48.tsp
cd ../java
./run_java.sh
```

## Testing
- Per-language `test/` subdirectories are for POC only - ignore them
- Future: Create centralized `tests/` directory for proper tests

## Performance Rule
- **No memory fragmentation**
- Use non-fragmented arrays/matrices as primary data structures

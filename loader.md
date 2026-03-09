# Loader

## Overview
The loader parses TSP (Traveling Salesman Problem) instance files and generates a distance matrix file required by the solvers.

## Location
`mlp-instances/loader/`

## Build
```bash
cd mlp-instances/loader
make
```

## Usage
You **must** be in the `mlp-instances/` directory to run the loader:

```bash
cd mlp-instances/
./load att48.tsp
```

This creates `distance_matrix` in `mlp/` (parent directory).

## Output Format
The `distance_matrix` file contains:

1. **Line 1**: Dimension (number of cities)
2. **Lines 2-N**: Upper-triangular distance matrix
3. **EOF**: End of matrix marker
4. **Filename**: Corresponding `.rnd` file name (e.g., `att48.rnd`)
5. **RND**: Random data marker
6. **Seed size**: Number of random values
7. **Random values**: Data from `.rnd` file

Example (showing end of file):
```
EOF
att48.rnd
RND
15745
12
3
1
...
```

## Dependencies
The loader requires `.rnd` files in `mlp-instances/rnd/`:
- Each `.tsp` instance needs a corresponding `.rnd` file
- Contains random seed data for GILS-RVND algorithm reproducibility
- If missing, loader aborts with: `Aborted: Ainda nao exieste arquivo '.rnd' para essa instancia.`

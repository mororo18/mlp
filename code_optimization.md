# Code Optimization

## Principle
- Performance optimization is **iterative** and **data-driven**
- Always **measure first**, then optimize

## Workflow Steps

1. **Understand** - Read and comprehend the implementation

2. **Baseline** - Run benchmarks to establish current performance

3. **Profile** - Identify bottlenecks (hot paths, memory issues)

4. **Optimize** - Ask permission before applying improvements

5. **Verify** - Re-benchmark to confirm improvements

6. **Document** - Ask permission to store metrics and explain what caused the improvement

7. **Update Implementation Details** - Update the language's `IMPLEMENTATION.md` file with changes made

---

## Implementation Details Files

- Each language directory should have an `IMPLEMENTATION.md` file
- Should describe implementation-specific details (data structures, algorithms, optimizations applied)
- Must be updated after any changes

---

## Profiling Tools

### Callgrind (Valgrind) - Main tool
```bash
valgrind --tool=callgrind ./program
```

**Tool needed**: Create `parse_callgrind.py` to extract relevant metrics from Callgrind output.

### perf (Linux)
```bash
perf record ./program
perf report
```

### Language-Specific
- Rust: `cargo flamegraph`
- Python: `cProfile`, `line_profiler`
- etc.

---

## Benchmark Usage

- `run_bm.py` - Run benchmarks
- `manager_bm.py` - Manage benchmark sessions

---

## Common Optimizations Examples

- **Memory layout**: Use contiguous arrays/matrices (core rule - no memory fragmentation)
- Reduce allocations in hot loops
- Cache-friendly data access patterns
- Language-specific micro-optimizations

# Table/overhead microbenchmarks

Standalone Lua scripts used to investigate (a) how Lua/LuaJIT actually
store `seq`-like flattened arrays internally (array part vs hash part of
a table), whether alternative data-structure shapes would be faster for
the access patterns the GILS-RVND hot loop actually uses, and (b) other
non-obvious runtime overheads (JIT compilation gaps, stdlib calls that
don't get inlined). Not part of the solver — run directly with
`lua5.3`/`luajit`, no dependencies.

Results and conclusions are written up in `../IMPLEMENTATION.md`. Kept
here (not in `/tmp`) so the investigation is reproducible and part of the
repo history, not just a conversation transcript.

- **`bench_01_array_vs_hash_part.lua`** — memory footprint + fill-time
  comparison (dense vs the real triangular `subseq_fill` gap pattern vs
  shuffled insertion order). Confirms Lua tables have a genuine
  array part (contiguous) separate from the hash part, and that
  **insertion order**, not density, decides which one an integer-keyed
  table lands in.
- **`bench_02_locality_square_vs_packed.lua`** — sequential read-locality
  test, current "square" (gapped) index vs a real triangular "packed"
  (no-gap) index, with both sides' index arithmetic equalized (row
  offsets precomputed) to isolate memory-layout effects from
  index-computation cost. (An earlier, unequalized version of this test
  gave a misleading result — see `IMPLEMENTATION.md` for that pitfall.)
- **`bench_03_struct_shapes_random.lua`** — small-scale (DIM=300) random-
  access comparison across 5 structure shapes: flat square, flat packed,
  nested `t[i][j][k]`, "array of structs" `t[i][j]={T=,C=,W=}`, and
  struct-of-arrays. Random access order defeats hardware prefetch, so
  it's a harsher test than sequential scan.
- **`bench_04_struct_shapes_intense.lua`** — large-scale (DIM=1500,
  ~100+ MB, well past this machine's 6 MiB L3) factorial test crossing
  {interleaved, struct-of-arrays} x {square, packed} — the flat survivors
  from bench_03 — under both full-random and windowed/neighborhood access
  (the latter mimicking `search_reinsertion`'s actual access shape).
- **`bench_05_table_insert_vs_manual.lua`** — `table.insert`/`table.remove`
  vs hand-written equivalents (`t[#t+1]=x`, manual shift-remove). Follows
  up on a `luajit -jv` finding: these two stdlib calls trigger trace
  "stitch" events (not inlined by the JIT) at their real call sites in
  `main.lua` (`construction()` and RVND's `neighbd_list` handling).
- **`bench_06_gc_pressure.lua`** — quantifies total allocation volume of a
  real solve (not a synthetic case) by running `main.lua` itself
  with `collectgarbage("stop")`, so nothing gets freed and the heap
  growth equals total bytes allocated. Must run from this `lua/`
  directory (`main.lua` resolves its own `dofile`/`io.open` paths
  relative to cwd). Answers whether `neighbd_list`'s per-improvement
  table literal (or any other per-iteration allocation) is a meaningful
  GC pressure source.
- **`bench_07_global_vs_local_stdlib.lua`** — global `math.floor`/
  `math.huge` access (as used throughout `main.lua`, never localized) vs
  `local floor = math.floor` hoisted once. Follows the LuaJIT Numerical
  Computing Performance Guide's advice to cache stdlib functions in
  upvalues.
- **`bench_08_ffi_vs_table.lua`** — LuaJIT-only (`ffi`, not available in
  lua5.3): flat `double[]` FFI array vs our current Lua-table layout, same
  DIM=1500 scale and access patterns as bench_04, to check whether the
  literature's often-cited large FFI speedups hold up against a layout
  that's already been optimized (they don't, at least not dramatically —
  see `IMPLEMENTATION.md`).
- **`bench_09_ffi_conditional_pattern.lua`** — validates the actual
  integration pattern for adopting FFI without a parallel luajit-only
  file: a single factory function branches on `pcall(require, "ffi")`,
  returning either an `ffi.new("double[?]", n+1)` array or a plain table;
  every other function keeps using `arr[idx]`/`t[idx]` unchanged since
  `[]` indexing is syntactically identical for both. Runs correctly on
  both lua5.3 (falls back to table) and luajit (uses FFI).

Run with e.g.:
```bash
lua5.3 bench_03_struct_shapes_random.lua
luajit bench_04_struct_shapes_intense.lua
```

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
  `main.canc.lua` (`construction()` and RVND's `neighbd_list` handling).

Run with e.g.:
```bash
lua5.3 bench_03_struct_shapes_random.lua
luajit bench_04_struct_shapes_intense.lua
```

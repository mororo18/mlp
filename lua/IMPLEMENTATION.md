# Lua — Implementation Notes

Decision log for `mlp-main/lua`. Format per `code_optimization.md` step 7.
Only `main` is optimized (see root `CLAUDE.md`) — this file documents
investigation and decisions for this branch only.

## Baseline (2026-07-08/09)

`run_bm.py --lang lua luajit -I lua-bench-instances.txt -n 10` against
`burma14`/`att48`/`rat99`. COST matched the expected value in every
repetition (20315 / 209320 / 57986), confirming the deterministic `.rnd`
trajectory is stable after the earlier sort-bug fix.

| instância | lua5.3 | luajit | speedup |
|---|---|---|---|
| burma14 (14) | ~0.14s | ~0.012s | ~11x |
| att48 (48) | ~18.0s | ~0.66s | ~27x |
| rat99 (99) | ~343s | ~12.2s | ~28x |

luajit's advantage grows with instance size, consistent with the JIT
compiling the hot loops that dominate at larger sizes.

## Profiling (2026-07-08)

Tool: `luajit -jp=a` (see `code_optimization.md` "Lua" section for why
this was chosen over `profiler.lua`/`ProFi.lua`). Run against `rat99`.

Aggregate by function (`-jp=f`):

| função | % tempo |
|---|---|
| `search_reinsertion` | 60% |
| `subseq_load` | 19% |
| `search_swap` | 10% |
| `search_two_opt` | 10% |

Source-annotated (`-jp=a`) pinpoints the hot lines inside
`search_reinsertion` (around lines 282, 340, 402, 424) and `subseq_load`
(~111-116) — all doing the same shape of access: reading `seq[...]` via
the flat index arithmetic `3 * ((s_size) * (x-1) + (y-1)) + z` (the old
`to_1D` macro, hand-expanded when the `canc` preprocessor dependency was
removed) combined with a distance-matrix lookup `c[x][y]`. Matches the
pattern seen in an old (2021, pre-refactor) `ProFi.lua` sample report
kept in this directory — `search_reinsertion` dominated there too, so
this isn't an artifact of the `rat99` instance specifically.

## Table-layout investigation (2026-07-09)

Question: is `seq` (and similar flattened structures) actually stored
contiguously by Lua, or is "everything a hash table" costing us
unnecessary overhead — and would restructuring the layout help the hot
loop above? Scripts in `bench/`, run on lua5.3 5.3.6 and LuaJIT 2.1
(this machine: 4 cores, L1d 32KB/core, L2 256KB/core, L3 6MB shared).

The testing methodology used here (random access order, confound
control, cache-busting scale, cheap-then-expensive funnel) is written up
language-agnostically in the root `code_optimization.md` ("Data-Structure
Layout Testing Methodology") so it can be reapplied to other languages —
this section is the Lua-specific run of it, not the methodology itself.

**Every conclusion below was validated by running the corresponding
script in this environment, not assumed from documentation.** One
intermediate test (see "confound" below) initially gave a wrong answer
and had to be corrected before being trusted — kept as a cautionary note
here rather than silently fixed, since it's the reason the following
results are trusted.

### 1. Lua tables do have a real array part (`bench_01`)

Confirmed empirically: an integer-keyed Lua table has a genuine
contiguous array part, separate from the hash part — and **insertion
order**, not key density, decides which one a given key lands in.
Filling a table with ascending (even gapped) integer keys — which is
what `subseq_fill`'s `for i=1,DIM do for j=i,DIM do ...` loop already
does — gets substantially better array-part promotion than filling the
same key set in shuffled order:

| padrão (DIM=100, 15150 chaves reais / 30000 slots) | lua5.3 KB/chave | luajit KB/chave |
|---|---|---|
| denso (1..N, sem buracos) | 0.0171 | 0.0086 |
| triangular (padrão real do `subseq_fill`) | 0.0254 (1.5x) | 0.0171 (~2x, footprint absoluto ≈ denso) |
| triangular embaralhado | 0.0423 (2.5x) | 0.0257 (3x) |

**Practical takeaway**: the current code's fill order is already
array-part-friendly. Don't reorder those loops.

### 2. Packing away the triangular gaps: negligible once arithmetic cost is controlled (`bench_02`)

First attempt at this comparison (not kept in `bench/`, see git history
of this file's early drafts if needed) computed the packed layout's row
offset with a division-heavy formula recomputed on *every* access, and
found packed ~3x slower — which turned out to be almost entirely the
cost of that division, not memory layout. Redone with both layouts'
per-row offset precomputed once (as real code would do), the difference
mostly disappeared:

- lua5.3: square 0.294s vs packed 0.291s (sequential full scan, DIM=700)
- luajit: square 0.021s vs packed 0.018s

**Practical takeaway**: eliminating `seq`'s triangular gaps is not a
meaningful win by itself. Don't pursue it in isolation.

### 3. What actually matters: avoid nested tables / per-entry named structs (`bench_03`)

Small-scale (DIM=300) random-access comparison (random order defeats
hardware prefetch, a harsher test than sequential scan) across 5 shapes
holding the same logical T/C/W triple per `(i,j)`:

| formato | lua5.3 | luajit |
|---|---|---|
| A) flat 1D square (código atual) | 0.0906s | 0.0122s |
| B) flat 1D packed triangular | 0.0774s (-15%) | 0.0128s (≈) |
| C) tabelas aninhadas `t[i][j][k]` | 0.1072s (+18%) | 0.0278s (**+128%**) |
| D) array of structs `t[i][j]={T=,C=,W=}` | 0.1190s (+31%) | 0.0386s (**+216%**) |
| E) struct of arrays (3 arrays flat) | 0.0750s (-17%) | 0.0124s (≈) |

**Practical takeaway**: the current code already avoids the expensive
shapes (C and D) by having flattened `seq` via the old `to_1D` macro
expansion. That flattening — not the specific square-vs-packed choice —
is where the real, validated win already lives. Confirms the project's
"no memory fragmentation, use contiguous arrays" rule empirically, not
just as a stated convention.

### 4. At realistic scale, the remaining flat-layout choices barely matter (`bench_04`)

Follow-up at DIM=1500 (~100+ MB per structure, well past this machine's
6 MiB L3, so any real cache effect should surface), factorial over the
flat survivors from step 3 ({interleaved, SoA} x {square, packed}),
under both full-random and a windowed/neighborhood access pattern meant
to mimic `search_reinsertion`'s actual bounded-range scanning:

**luajit** (the interpreter that matters for `main`, ~28x faster than
lua5.3 at `rat99` scale):
- Full random: interleaved+square (current) fastest at 0.42s; the other
  three cluster at 0.46-0.48s (~10-14% slower).
- Windowed: SoA+square fastest at 0.0094s; interleaved+square close
  second at 0.0111s; the other two within ~35% of the fastest.

lua5.3 showed a larger relative swing (interleaved winning full-random
by ~28%, SoA winning windowed by ~22%), but luajit — the interpreter
that actually matters here — keeps the current layout competitive or
best in both patterns tested.

**Conclusion**: no restructuring of `seq`'s layout (packing, SoA) is
clearly justified by this data for luajit. The current flat/interleaved
layout is already a reasonable choice. Not pursuing further changes to
`seq`'s shape based on this investigation.

## Open follow-up (not started)

Registered in root `TODO.md`: investigate other non-obvious overhead
sources in this implementation beyond table layout (function-call/
closure overhead, `table.remove`/`table.insert` cost inside RVND,
`--` string comparisons, GC pauses, etc.) — this investigation only
looked at `seq`'s data-structure shape, not the broader hot-loop cost
breakdown implied by the 19%+60% profiler split above.

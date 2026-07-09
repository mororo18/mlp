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

## Why is lua5.3 so much slower than luajit? (2026-07-09)

Decomposed empirically rather than assumed. LuaJIT has a `-joff` flag
that disables the JIT compiler entirely, leaving only its bytecode
interpreter running — comparing that against normal luajit isolates "the
interpreter is better" from "the JIT compiler helps", and comparing it
against lua5.3 isolates interpreter-vs-interpreter architecture from
compilation. Same `att48` instance, three configurations:

```bash
lua5.3 main.canc.lua        # 18.54s -- PUC-Lua's C-switch interpreter
luajit -joff main.canc.lua  #  6.37s -- LuaJIT's interpreter, JIT disabled
luajit main.canc.lua        #  0.53s -- LuaJIT's interpreter + JIT compiler
```

Two independent, multiplicative factors, not one:

1. **~2.9x from the interpreter alone** (18.54s → 6.37s), with *zero*
   compilation happening. LuaJIT's interpreter is hand-written in
   assembly (per-architecture, via DynAsm) with a tight dispatch loop,
   vs lua5.3's portable-C switch-based interpreter — plus LuaJIT uses
   NaN-tagged 8-byte values (numbers/bools/pointers unboxed into the
   unused bits of a NaN double) vs Lua 5.3's larger tagged-union
   `TValue`, meaning less memory traffic per operation even before any
   JIT is involved. External corroboration: "LuaJIT's hand-coded-in-
   assembly interpreter... running almost 3x the speed" of PUC-Lua's —
   matches the measured 2.9x almost exactly. ([Building the fastest Lua
   interpreter, automatically!](https://sillycross.github.io/2022/11/22/2022-11-22/))
2. **~12x more from the JIT compiler** (6.37s → 0.53s) — the tracing
   compiler identifies the hot loop (`search_reinsertion` etc., see
   Profiling above), specializes it for the types actually observed, and
   emits native machine code, skipping bytecode dispatch entirely for
   those loops.

Total: ~35x (18.54s / 0.53s), decomposed as 2.9x × 12x ≈ 35x. Neither
factor alone explains the gap — an interpreter-only comparison
(lua5.3 vs `luajit -joff`) would undersell the JIT's contribution, and
citing "LuaJIT has a JIT compiler" alone ignores that its interpreter is
already ~3x faster before compilation ever kicks in.

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

## Non-obvious overhead investigation (2026-07-09)

Follow-up to the table-layout work above: what else, beyond `seq`'s data-
structure shape, might be costing real time in this implementation
without showing up as an obvious line in the `-jp=a` output? Method: run
LuaJIT's own compilation log (`luajit -jv`) against the real
`main.canc.lua` solving `rat99`, which shows every trace compiled, every
trace *abort* (falls back to slow interpretation entirely) and every
trace *stitch* (bails out to run one uncompiled operation via the
interpreter, then resumes — LuaJIT 2.1's mitigation for NYI/"not yet
implemented" stdlib functions, cheaper than an abort but not free).

```bash
luajit -jv main.canc.lua > /dev/null 2> jv.log
grep -c abort jv.log     # 0
grep -c stitch jv.log    # 22
grep stitch jv.log | grep -oE "stitch [a-zA-Z_.]+" | sort | uniq -c
#   7 stitch table.insert
#   4 stitch string.gmatch_aux
#   1 stitch string.gmatch
```

### 1. Zero trace aborts — no hidden "falls back to slow interpreter" disaster

Good news, worth stating explicitly rather than assuming: every hot loop
(`search_reinsertion`, `subseq_load`, `search_swap`, `search_two_opt`)
compiles to native code with no full abort. Nothing here is silently
running fully interpreted.

### 2. `table.insert`/`table.remove` are not JIT-inlined — confirmed at their real call sites

7 of the 22 stitch events are `table.insert`; these trace back to
`construction()` (lines 186-220) and RVND's `neighbd_list` handling
(line 512, `table.remove(neighbd_list, index)`). Corroborated externally:
LuaJIT 2.1 added trace stitching specifically to handle NYI stdlib calls
like these without a full abort ([LuaJIT NYI
list](https://github.com/tarantool/tarantool/wiki/LuaJIT-Not-Yet-Implemented),
[Percona: The Anatomy of LuaJIT
Tables](https://percona.community/blog/2020/04/29/the-anatomy-of-luajit-tables-and-whats-special-about-them/)) —
expected LuaJIT behavior, not a peculiarity of this environment.

Validated the fix empirically (`bench_05_table_insert_vs_manual.lua`):
hand-written equivalents are consistently faster in both interpreters,
though the margin is modest, not dramatic:

| operação | lua5.3 | luajit |
|---|---|---|
| `table.insert(t,x)` | 0.362s | 0.0138s |
| manual `t[#t+1]=x` | 0.283s (-22%) | 0.0124s (-10%) |
| `table.remove(list,3)` | 0.138s | 0.0475s |
| manual shift-remove | 0.100s (-27%) | 0.0440s (-7%) |

**Practical takeaway**: replacing `table.insert`/`table.remove` with
hand-written equivalents in `construction()` and RVND's `neighbd_list`
handling is a validated, low-risk, modest-magnitude optimization
candidate — this is Optimize-phase work, not yet applied, needs
permission per the workflow in `code_optimization.md`.

### 3. `string.gmatch` stitches are confined to one-time file parsing — not a hot-loop concern

The other 5 stitch events (`string.gmatch`/`string.gmatch_aux`) trace to
`Data.lua`'s instance-file parsing, which runs once at startup. Per the
general LuaJIT guidance found alongside the NYI docs above: avoiding NYI
calls only matters for code that runs enough times for the JIT to want to
compile it — a one-shot initialization routine is never a target for
compilation regardless, so this is correctly not a concern.

### 4. No closures in the hot path

Checked directly (`grep "function("` across `main.canc.lua`): the only
anonymous function in the file is inside `protect()`, which is dead code
(its one call site is commented out, line 704). No closure-allocation
overhead exists in the current hot loop.

### 5. GC pressure — investigated and ruled out (`bench_06`, 2026-07-09)

`neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}` (line
509) allocates a fresh 5-element table every time RVND finds an improving
move — flagged as a plausible GC churn source, now quantified. Method:
run the real `main.canc.lua` solve with `collectgarbage("stop")`
beforehand, so nothing gets freed — the heap growth then equals total
bytes *allocated* (garbage + live), not just what's currently live.

| instância | interpretador | tempo | total alocado (GC desligado) |
|---|---|---|---|
| att48 (48) | lua5.3 | 18.0s | 2.7 MB |
| att48 (48) | luajit | 0.50s | 2.0 MB |
| rat99 (99) | luajit | 13.0s | 8.2 MB |

**Conclusion**: ruled out. Going from `att48` to `rat99`, allocation
volume grows ~4x (2.0→8.2 MB) while runtime grows ~26x (0.5s→13.0s) —
allocation does not scale anywhere close to proportionally with runtime,
so it cannot be what dominates the cost at larger instance sizes. A few
MB of total garbage over a multi-second run is trivial for any GC to
collect (low milliseconds at most) — nowhere near enough to explain a
meaningful fraction of the measured time. `neighbd_list` and friends are
not a genuine optimization target. This closes the item registered in
root `TODO.md`.

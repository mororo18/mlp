# TODO

## Usage

Every task here should be criticized and not taken as gospel. They may be out of date or misguided in the first place. You should thoroughly consider them, do research, possibly push back, and ask clarifying questions.

**Priorities are section-based:** P1 highest, then P2, P3, and so on.

Sometimes within a section, tasks may be grouped (extra newline separating them from other groups) to imply they should be considered together.

**Current task:** first item or group in the highest-priority non-empty section.

**Multi-language tasks:** When a task applies to every language implementation, track which languages have been completed directly in the TODO file. Apply the tracking immediately as each language is completed, not after finishing all languages.

**Planning:** Always make a plan before trying to do a task and ask for approval before executing it.

**Findings:** Add relevant findings discovered when working on a task to the TODO file.

**Completion:** When a task is finished, remove it from the TODO file rather than marking it as done.

**Commits:** Do not automatically commit and push changes after completing a task. Ask for approval before committing.

---

## P2

Fix `lua/run_lua.sh` and `lua/run_luajit.sh`: they invoke `main.canc.lua.lua`,
which crashes under both `lua5.3` (`table.remove` out of bounds in `RVND`,
main.canc.lua.lua:343) and `luajit` (`attempt to index a nil value` in
`sort`, main.canc.lua.lua:91) — reproduced on a clean checkout, unrelated to
any pending change. Several near-duplicate files in the same directory
(`main.lua`, `main_lua.lua`, `main_luajit.canc`, `main_luajit.canc.lua`) run
correctly to completion under `luajit`, but weren't verified under `lua5.3`.
Needs a decision on which file is canonical (ties into the item below) before
`run_lua.sh`/`run_luajit.sh` can be fixed to point at a working entry point.

Blocks: the lua leg of the cross-branch progress-flag rollout
(`tcc/TODO.md` P2) — skipped for now pending this fix.

Fix `octave/main.m` and `octave/Data.m`: never actually runs to completion
on a clean checkout, two independent pre-existing bugs found while trying
to verify it before adding the progress flag:
1. `main.m` calls `mainn()` on line 2, before the `mainn` function is
   defined further down in the same script file — Octave errors
   `'mainn' undefined` immediately. Local functions in an Octave script
   must be defined before they're invoked at the top level; the call
   needs to move to the end of the file, after all function definitions.
2. Once that's fixed, `octave/Data.m:9` does
   `str2double(lines{index}(:))` — the `(:)` reshapes the dimension
   string (e.g. `"48"`) into a column vector of individual characters,
   so `str2double` parses each digit separately (`[4;8]`) instead of the
   whole number, corrupting `dimension` (`NaN`) for any instance with a
   2+ digit dimension — i.e. any real instance. Needs a proper
   string-to-number parse (e.g. plain `str2double(lines{index})` without
   the `(:)`).
   Also worth a look: `main.m`'s local `sort` function shadows Octave's
   builtin `sort`, which broke a builtin-internal call
   (`ismember`/`​__unimplemented__`) once `mainn` was reachable — renaming
   the local helper (e.g. to `local_sort`) avoids the collision.

Blocks: the octave leg of the cross-branch progress-flag rollout
(`tcc/TODO.md` P2) — skipped for now pending this fix. `run_matlab.sh`
uses the same `main.m`, so is presumably affected too, but couldn't be
verified (`matlab` unavailable in this environment).

Fix `csharp` mono build path: `mcs -debug -optimize+ -out:solve_mcs *.cs`
(used by `build.sh`/`run_mcs.sh`) fails with
`GILS_RVND.cs(87,13): error CS0131` — this Mono compiler (6.12.0.199)
doesn't support the C# 7 tuple-deconstruction assignment syntax
(`(arr[i + 1], arr[right]) = (arr[right], arr[i + 1]);`) used in
`Partition`. Reproduced on a clean checkout, unrelated to any pending
change. The `dotnet build` path works fine and was used to verify the
progress-flag addition. Separately, running the built `net6.0` binary
directly (`./bin/Debug/net6.0/csharp`, what `run_dotnet.sh` calls) fails
in this environment with a missing-framework error since only the 8.0
runtime is installed side-by-side; `DOTNET_ROLL_FORWARD=LatestMajor`
works around it here, but the underlying net6.0-target-vs-8.0-installed
mismatch (see `SETUP.md`'s Tool Versions table) may need a decision
(retarget the csproj to net8.0, or install net6.0 runtime).

Review multiple main variants in `lua/`, `python/` - keep canonical version

Remove profiler output files (`.txt`, `.log`, `MyProfilingReport.txt`)

Remove compiled binaries from repo

---

## P4

Create `benchmarks/` directory for storing metrics

Create `parse_callgrind.py` tool to extract relevant metrics from Callgrind output

Create proof of concept small benchmark (if not existent) that shows the change is really relevant

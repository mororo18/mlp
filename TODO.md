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

Review multiple main variants in `lua/`, `python/` - keep canonical version

Remove profiler output files (`.txt`, `.log`, `MyProfilingReport.txt`)

Remove compiled binaries from repo

---

## P4

Create `benchmarks/` directory for storing metrics

Create `parse_callgrind.py` tool to extract relevant metrics from Callgrind output

Create proof of concept small benchmark (if not existent) that shows the change is really relevant

# TODO

## Usage

Every task here should be criticized and not taken as gospel. They may be out of date or misguided in the first place. You should thoroughly consider them, do research, possibly push back, and ask clarifying questions.

**Priorities are section-based:** P1 highest, then P2, P3, and so on.

Sometimes within a section, tasks may be grouped (extra newline separating them from other groups) to imply they should be considered together.

**Current task:** first item or group in the highest-priority non-empty section.

---

## P1

Remove local `loader/` directory (deprecated, relevant loader is in `mlp-instances/loader/`)

---

## P2

Create centralized `tests/` directory

Remove per-language `test/` subdirectories (POC tests)

---

## P3

Ensure every language dir has `build.sh` and `run_*.sh`

---

## P4

Review multiple main variants in `lua/`, `python/` - keep canonical version

Remove profiler output files (`.txt`, `.log`, `MyProfilingReport.txt`)

Remove compiled binaries from repo

---

## P5

Create `benchmarks/` directory for storing metrics

Create `parse_callgrind.py` tool to extract relevant metrics from Callgrind output

Create proof of concept small benchmark (if not existent) that shows the change is really relevant

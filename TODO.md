# TODO

## Usage

Every task here should be criticized and not taken as gospel. They may be out of date or misguided in the first place. You should thoroughly consider them, do research, possibly push back, and ask clarifying questions.

**Priorities are section-based:** P1 highest, then P2, P3, and so on.

Sometimes within a section, tasks may be grouped (extra newline separating them from other groups) to imply they should be considered together.

**Current task:** first item or group in the highest-priority non-empty section.

**Multi-language tasks:** When a task applies to every language implementation, track which languages have been completed directly in the TODO file. Apply the tracking immediately as each language is completed, not after finishing all languages.

**Planning:** Always make a plan before trying to do a task and ask for approval before executing it.

**Completion:** When a task is finished, remove it from the TODO file rather than marking it as done.

---

## P2

Create centralized `tests/` directory

Remove per-language `test/` subdirectories (POC tests)

---

## P3

Ensure every language dir has `build.sh` and `run_*.sh`

Add progress output flag to all language implementations (default: false)

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

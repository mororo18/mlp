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

Add progress output flag to all language implementations (default: false)

**Note:** Before starting work on each language, verify the implementation builds and runs successfully.

- java [DONE]
- rust [DONE]
- python [DONE]
- c [DONE]
- cplusplus
- python
- c
- cplusplus
- fortran
- go
- javascript
- julia
- lua
- octave
- csharp
- c_asm
- cppOOP

---

## P3

Review multiple main variants in `lua/`, `python/` - keep canonical version

Remove profiler output files (`.txt`, `.log`, `MyProfilingReport.txt`)

Remove compiled binaries from repo

---

## P4

Create `benchmarks/` directory for storing metrics

Create `parse_callgrind.py` tool to extract relevant metrics from Callgrind output

Create proof of concept small benchmark (if not existent) that shows the change is really relevant

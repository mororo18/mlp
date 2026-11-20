#!/usr/bin/env python3
"""Parallel, resumable benchmark runner.

Shared machinery (sandboxing, distance-matrix cache, resume, memory
sampling) lives in core.py -- see its module docstring for the race
condition that makes the sandboxing necessary. This script benchmarks
whatever each language's build.sh already produces. For a debug-vs-release
comparison, see bench_modes.py.

-t/--translations also benchmarks LLM translations (one directory per
<from>_<to> project). Their build and run commands come from the PROJETOS
catalog in <translations dir>/benchmark.py, whose toolchains replicate each
reference's build.sh flags. Results go to <project>.csv, `source` = project
name (cpp_lua's LuaJIT run is "cpp_lua_luajit").
"""

import argparse
import concurrent.futures
import dataclasses
import importlib.util
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import core

REPO_ROOT = core.REPO_ROOT
DEFAULT_OUT = REPO_ROOT.parent / "mlp_result"
SOURCES = core.SOURCES
LANG_DIR = core.LANG_DIR
RUN_SCRIPT = core.RUN_SCRIPT
DEFAULT_TRANSLATIONS_DIR = Path.home() / "Downloads" / "analise" / "claude"


@dataclasses.dataclass(frozen=True)
class Target:
    name: str          # `source` column
    csv: str           # output CSV name, also the sandbox subdirectory
    src_dir: Path
    run: tuple = None  # None: parse the language's run_*.sh


def lang_target(lang):
    return Target(lang, LANG_DIR[lang], REPO_ROOT / LANG_DIR[lang])


# --------------------------------------------------------------------------
# One-time setup: per-language build.sh
# --------------------------------------------------------------------------

def build_languages(langs):
    done = set()
    for lang in langs:
        ld = LANG_DIR[lang]
        if ld in done:
            continue
        done.add(ld)
        build_script = REPO_ROOT / ld / "build.sh"
        if build_script.exists():
            core.log(f"[build] {ld}")
            subprocess.run(["bash", "build.sh"], cwd=REPO_ROOT / ld, check=True)


def load_translation_catalog(tdir):
    spec = importlib.util.spec_from_file_location("translations_benchmark", tdir / "benchmark.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.ENV = mod.carregar_ambiente()  # its builders run with this env
    return mod


def translation_targets(mod, tdir, projects):
    targets, builds = [], []
    for proj in projects:
        var = mod.PROJETOS[proj]["traduzida"]
        targets.append(Target(proj, proj, tdir / proj, tuple(var["run"])))
        for extra in var.get("runtimes_extras", []):
            targets.append(Target(f"{proj}_{extra['lang'].lower()}", proj, tdir / proj,
                                  tuple(extra["run"])))
        if var["build"]:
            builds.append((proj, tdir / proj, var["build"]))
    return targets, builds


def build_translations(mod, builds):
    for proj, cwd, (toolchain, spec) in builds:
        core.log(f"[build] {proj}")
        ok, err, _ = mod.TOOLCHAINS[toolchain](str(cwd), spec)
        if not ok:
            sys.exit(f"[build] {proj} failed:\n{err}")


def get_run_command(target, sandbox_lang_dir):
    if target.run:
        return list(target.run)
    return core.run_command_from_script(sandbox_lang_dir / RUN_SCRIPT[target.name])


# --------------------------------------------------------------------------
# Task scheduling
# --------------------------------------------------------------------------

@dataclasses.dataclass
class Task:
    target: Target
    instance: str
    rep: int


def build_tasks(targets, instances, n, counts, fresh):
    tasks = []
    for inst in instances:
        for t in targets:
            done = 0 if fresh else counts.get((t.name, inst), 0)
            for rep in range(max(0, n - done)):
                tasks.append(Task(t, inst, rep))
    return tasks


timed_out_langs = set()
timed_out_lock = threading.Lock()


def process_task(task, timeout, sandbox_root, dm_cache, writer, branch, cores):
    name = task.target.name
    with timed_out_lock:
        if name in timed_out_langs:
            core.log(f"[skip] {name} already timed out this run, skipping {task.instance}")
            return

    slot = core.get_slot_id()
    sandbox = core.Sandbox(sandbox_root, slot)
    lang_path = sandbox.lang_dir(task.target.csv, task.target.src_dir)
    sandbox.set_distance_matrix(dm_cache[task.instance])

    cpu_core = cores[slot % len(cores)] if cores else None

    core.log(f"[slot{slot}] start {name:<8} {task.instance:<16} rep={task.rep}"
             + (f" core={cpu_core}" if cpu_core is not None else ""))
    t0 = time.time()
    cmd = get_run_command(task.target, lang_path)
    info = core.run_once(cmd, lang_path, timeout, cpu_core)
    dt = time.time() - t0
    info.update({"source": name, "instance": task.instance, "branch": branch})
    writer.append(task.target.csv, info)

    if info["timed_out"]:
        with timed_out_lock:
            timed_out_langs.add(name)
        core.log(f"[slot{slot}] TIMEOUT {name} {task.instance} after {timeout}s "
                 "- skipping its remaining queued runs")
    else:
        core.log(f"[slot{slot}] done  {name:<8} {task.instance:<16} rep={task.rep} "
                 f"COST={info['COST']} TIME={info['TIME']} ({dt:.1f}s wall)")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_arg_parser():
    p = argparse.ArgumentParser(
        description="Run MLP benchmarks across languages, in parallel, with resume support.",
        parents=[core.common_parser(DEFAULT_OUT)])
    p.add_argument("-l", "--lang", nargs="+",
                    help=f"Languages to test, or 'all' (default when -t is not given). "
                         f"Choices: {', '.join(SOURCES)}")
    p.add_argument("-t", "--translations", nargs="+",
                    help="Translation projects to test (e.g. cpp_go java_cpp), or 'all'")
    p.add_argument("--translations-dir", default=str(DEFAULT_TRANSLATIONS_DIR),
                    help="Directory with the <from>_<to> projects and benchmark.py")
    return p


def resolve_langs(requested):
    if requested is None:
        return []
    if requested == ["all"]:
        return SOURCES[:]
    for lang in requested:
        if lang not in SOURCES:
            sys.exit(f"{lang} is not supported. Choose from: {', '.join(SOURCES)}")
    return requested


def main():
    args = build_arg_parser().parse_args()

    langs = resolve_langs(args.lang if args.lang or args.translations else ["all"])
    targets = [lang_target(l) for l in langs]
    builds = []
    if args.translations:
        tdir = Path(args.translations_dir).resolve()
        mod = load_translation_catalog(tdir)
        projects = sorted(mod.PROJETOS) if args.translations == ["all"] else args.translations
        unknown = [p for p in projects if p not in mod.PROJETOS]
        if unknown:
            sys.exit(f"Unknown translation(s): {', '.join(unknown)}. "
                     f"Choose from: {', '.join(sorted(mod.PROJETOS))}")
        extra, builds = translation_targets(mod, tdir, projects)
        targets += extra
    instances = core.resolve_instances(args)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    branch = core.get_branch()

    core.build_loader()
    dm_cache = core.precompute_distance_matrices(instances, core.SCRATCH / "dm_cache")
    build_languages(langs)
    if builds:
        build_translations(mod, builds)

    counts = core.existing_counts(out_dir, sorted({t.csv for t in targets}))
    tasks = build_tasks(targets, instances, args.repeat, counts, args.fresh)

    if not tasks:
        core.log("Nothing to do: all requested runs are already present in --out "
                  "(use --fresh to force a rerun).")
        return

    core.log(f"{len(tasks)} run(s) queued across {len(targets)} target(s), "
             f"{len(instances)} instance(s), {args.jobs} parallel job(s)")

    sandbox_root = core.SCRATCH / "sandbox"
    if sandbox_root.exists():
        shutil.rmtree(sandbox_root)
    sandbox_root.mkdir(parents=True)

    cores = core.resolve_pin_and_cores(args)

    writer = core.ResultWriter(out_dir)

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [
                pool.submit(process_task, t, args.timeout, sandbox_root, dm_cache, writer, branch, cores)
                for t in tasks
            ]
            for f in concurrent.futures.as_completed(futures):
                f.result()
    except KeyboardInterrupt:
        core.log("\n[interrupted] results recorded so far are saved; rerun the same "
                 "command to resume.")
        raise
    finally:
        if not args.keep_sandbox:
            shutil.rmtree(sandbox_root, ignore_errors=True)

    core.log("Benchmark run complete.")


if __name__ == "__main__":
    main()

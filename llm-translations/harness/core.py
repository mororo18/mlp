"""Shared machinery for bench.py and bench_modes.py.

Every solver reads its instance data from a hardcoded relative path,
"../distance_matrix" (see e.g. c/data.c, java/Data.java, python/read.py),
and mlp-instances/load always (re)writes that single file at the repo
root. That makes naive parallel execution unsafe: two runs against
different instances at the same time would race on that one file.

To parallelize safely without touching any language's source, each
concurrent "slot" gets a private sandbox directory that mirrors a
language directory (via symlinks to the real files/binaries) sitting
next to its own private distance_matrix. Solvers run with that sandbox
as their cwd, so "../distance_matrix" resolves to a per-slot file
instead of the shared one.
"""

import argparse
import csv
import itertools
import os
import shlex
import shutil
import subprocess
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd
import psutil

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRATCH = REPO_ROOT / ".bench"
DEFAULT_INSTANCES_FILE = REPO_ROOT / "inst-list"
TEMPLATE_CSV = REPO_ROOT / "data" / "template.csv"
DEFAULT_TIMEOUT = 60 * 60 * 24 * 2  # seconds per individual run

LANG_DIR = {
    "c": "c",
    "cplusplus": "cplusplus",
    "fortran": "fortran",
    "go": "go",
    "rust": "rust",
    "swift": "swift",
    "julia": "julia",
    "java": "java",
    "csharp": "csharp",
    "javascript": "javascript",
    "pypy": "python",
    "luajit": "lua",
    "lua": "lua",
    "python": "python",
}

RUN_SCRIPT = {
    "c": "run_c.sh",
    "cplusplus": "run_cpp.sh",
    "fortran": "run_fortran.sh",
    "go": "run_golang.sh",
    "rust": "run_rust.sh",
    "swift": "run_swift.sh",
    "julia": "run_julia.sh",
    "java": "run_java.sh",
    "csharp": "run_dotnet.sh",
    "javascript": "run_node.sh",
    "pypy": "run_pypy.sh",
    "luajit": "run_luajit.sh",
    "lua": "run_lua.sh",
    "python": "run_python3.sh",
}

SOURCES = list(LANG_DIR)

print_lock = threading.Lock()


def log(msg):
    with print_lock:
        print(msg, flush=True)


def get_branch():
    out = subprocess.check_output(["git", "branch"], cwd=REPO_ROOT).decode()
    for line in out.splitlines():
        if line.startswith("*"):
            return line[1:].strip()
    return ""


# --------------------------------------------------------------------------
# Result CSVs
# --------------------------------------------------------------------------

class ResultWriter:
    """Thread-safe, append-only writer, one lock per output CSV file."""

    def __init__(self, out_dir, template=TEMPLATE_CSV):
        self.out_dir = out_dir
        self.template = template
        self._locks = {}
        self._meta_lock = threading.Lock()

    def _lock_for(self, csv_name):
        with self._meta_lock:
            if csv_name not in self._locks:
                self._locks[csv_name] = threading.Lock()
            return self._locks[csv_name]

    def append(self, csv_name, row):
        path = self.out_dir / f"{csv_name}.csv"
        with self._lock_for(csv_name):
            if not path.exists():
                shutil.copy(self.template, path)
            with open(path) as f:
                header = f.readline().strip().split(",")
            clean_row = {k: ("" if v is None else v) for k, v in row.items()}
            with open(path, "a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=header, extrasaction="ignore", restval="")
                writer.writerow(clean_row)


def existing_counts(out_dir, csv_names, key_cols=("source", "instance")):
    """key_cols tuple -> number of rows already recorded, for resume."""
    key_cols = list(key_cols)
    counts = {}
    for name in csv_names:
        path = out_dir / f"{name}.csv"
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            continue
        if not all(c in df.columns for c in key_cols):
            continue
        for key, n in df.groupby(key_cols).size().items():
            key = key if isinstance(key, tuple) else (key,)
            counts[key] = counts.get(key, 0) + int(n)
    return counts


# --------------------------------------------------------------------------
# One-time setup: loader build, distance_matrix cache
# --------------------------------------------------------------------------

def needs_build(build_dir, rebuild):
    """True if the caller must (re)build into build_dir: forced via
    `rebuild`, the directory doesn't exist yet, or a previous attempt there
    never finished (no .built marker -- e.g. it crashed partway through a
    copytree/build, which used to be silently treated as "already built" on
    the next run). Call mark_built() right after a successful build."""
    marker = build_dir / ".built"
    if build_dir.exists():
        if not rebuild and marker.exists():
            return False
        shutil.rmtree(build_dir)
    return True


def mark_built(build_dir):
    (build_dir / ".built").touch()


def build_loader():
    log("[setup] building mlp-instances/loader")
    subprocess.run(["make"], cwd=REPO_ROOT / "mlp-instances" / "loader", check=True)


def precompute_distance_matrices(instances, cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    mlp_instances_dir = REPO_ROOT / "mlp-instances"
    mapping = {}
    for inst in instances:
        cache_path = cache_dir / (inst.replace("/", "_") + ".distance_matrix")
        if not cache_path.exists():
            log(f"[setup] generating distance matrix for {inst}")
            subprocess.run(["./load", inst], cwd=mlp_instances_dir, check=True)
            shutil.copy(REPO_ROOT / "distance_matrix", cache_path)
        mapping[inst] = cache_path
    return mapping


# --------------------------------------------------------------------------
# Per-slot sandbox
# --------------------------------------------------------------------------

class Sandbox:
    """Owned by exactly one worker thread for its whole lifetime (see
    get_slot_id), so no locking is needed for its internal state."""

    def __init__(self, root, slot):
        self.root = root / f"slot{slot}"
        self.root.mkdir(parents=True, exist_ok=True)
        self._populated = set()

    def lang_dir(self, name, source_dir):
        """Mirror source_dir's top-level entries as symlinks under a real
        (non-symlink) sandbox directory called `name`. Must stay a real
        directory: solvers cd into it and read "../distance_matrix", and
        that ".." resolution happens against the physical path, not any
        logical symlink chain."""
        target = self.root / name
        if name not in self._populated:
            target.mkdir(parents=True, exist_ok=True)
            for entry in os.scandir(source_dir):
                link = target / entry.name
                if not link.exists() and not link.is_symlink():
                    os.symlink(Path(entry.path).resolve(), link)
            self._populated.add(name)
        return target

    def set_distance_matrix(self, cache_path):
        link = self.root / "distance_matrix"
        if link.exists() or link.is_symlink():
            link.unlink()
        os.symlink(cache_path.resolve(), link)


# --------------------------------------------------------------------------
# Running a single benchmark
# --------------------------------------------------------------------------

def get_mem_avg(pid, timeout):
    p = psutil.Process(pid=pid)
    start = time.time()
    MB = 1024 ** 2
    gap = 1e-10

    mem_max = -99999
    mem_avg = 0
    count = 0
    timed_out = False
    while True:
        if time.time() - start >= timeout:
            timed_out = True
            try:
                p.kill()
            except psutil.NoSuchProcess:
                pass
            break

        try:
            with p.oneshot():
                if p.status() == psutil.STATUS_ZOMBIE:
                    break
                mem = p.memory_info()
                mem_total = (float(mem.rss) - float(mem.shared)) / MB
                if mem_total < gap:
                    break
                mem_max = max(mem_max, mem_total)
                mem_avg += mem_total
        except psutil.NoSuchProcess:
            break
        count += 1
        time.sleep(0.05)

    return {
        "mem_avg": mem_avg / count if count else float("nan"),
        "mem_max": mem_max if count else float("nan"),
        "mem_lookups": count,
        "timed_out": timed_out,
    }


def parse_metric(output, key):
    for line in output.splitlines():
        if key in line:
            i = line.find(":")
            return float(line[i + 1:].replace(",", "."))
    return None


def run_command_from_script(script_path):
    """Every run_*.sh has a 2-line shape: shebang, then the actual
    invocation. Measuring memory needs the solver's own pid, not a bash
    wrapper's, so the command is parsed out and run directly."""
    with open(script_path) as f:
        f.readline()  # shebang
        line = f.readline().strip()
    return [tok for tok in shlex.split(line) if tok != "$@"]


def run_once(cmd, cwd, timeout, cpu_core=None):
    with tempfile.TemporaryFile() as tempf:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=tempf, stderr=subprocess.STDOUT)
        if cpu_core is not None:
            try:
                psutil.Process(proc.pid).cpu_affinity([cpu_core])
            except psutil.Error:
                pass

        info = get_mem_avg(proc.pid, timeout)

        proc.wait()
        tempf.seek(0)
        output = tempf.read().decode("utf-8", errors="replace")

    if info["timed_out"]:
        info["COST"], info["TIME"] = None, None
    else:
        info["COST"] = parse_metric(output, "COST")
        wall = parse_metric(output, "wall clock")
        info["TIME"] = wall if wall is not None else parse_metric(output, "TIME")
    return info


# --------------------------------------------------------------------------
# Slot assignment
# --------------------------------------------------------------------------

_thread_local = threading.local()
_slot_counter = itertools.count()


def get_slot_id():
    if not hasattr(_thread_local, "slot"):
        _thread_local.slot = next(_slot_counter)
    return _thread_local.slot


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def resolve_instances(args):
    if args.instance:
        return [args.instance]
    with open(args.instances_list) as f:
        return [line.strip() for line in f if line.strip()]


def common_parser(default_out, default_timeout=DEFAULT_TIMEOUT,
                   default_instances=DEFAULT_INSTANCES_FILE):
    """Flags shared by bench.py and bench_modes.py. Each script adds its
    own -l/--lang (choices differ) and any script-specific flags on top,
    via parents=[common_parser(...)]."""
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i", "--instance", help="Run a single instance")
    p.add_argument("-I", "--instances-list", default=str(default_instances),
                    help="File with one instance per line (default: inst-list)")
    p.add_argument("-n", "--repeat", type=int, default=1,
                    help="Repetitions per (language, instance)")
    p.add_argument("-j", "--jobs", type=int, default=1,
                    help="Number of benchmark runs to execute concurrently")
    p.add_argument("--timeout", type=int, default=default_timeout,
                    help="Per-run timeout in seconds")
    p.add_argument("--out", default=str(default_out), help="Output directory for result CSVs")
    p.add_argument("--fresh", action="store_true",
                    help="Ignore existing results and rerun everything (default: resume)")
    pin = p.add_mutually_exclusive_group()
    pin.add_argument("--pin", dest="pin", action="store_true", default=None,
                       help="Pin each parallel job to its own CPU core (default: on when --jobs > 1)")
    pin.add_argument("--no-pin", dest="pin", action="store_false",
                       help="Disable CPU pinning")
    p.add_argument("--keep-sandbox", action="store_true",
                    help="Keep per-job sandbox directories after the run (debugging)")
    return p


def resolve_pin_and_cores(args):
    """Resolve the --pin default and, if pinning, the list of usable cores."""
    if args.pin is None:
        args.pin = args.jobs > 1
    cores = None
    if args.pin:
        cores = psutil.Process().cpu_affinity()
        if len(cores) < args.jobs:
            log(f"[warn] only {len(cores)} usable core(s) for {args.jobs} job(s); "
                "some jobs will share a core")
    return cores

"""Monta as variantes do experimento e o catálogo que o bench.py -t lê.

Cada linha de variants/manifest.csv é uma variante: uma cópia de
translations/<base>/ (ou da nossa implementação, com base ours/<pasta>) com os patches da coluna `patches` aplicados em ordem
(caminhos relativos a variants/, separados por espaço). Nada é editado à mão.

Saída em --out: uma pasta por variante e um benchmark.py com PROJETOS,
TOOLCHAINS e carregar_ambiente(), que é a interface que o bench.py espera.
A configuração de build e execução de cada variante é a da sua base no
catálogo do harness antigo (harness/catalog_base.py), com as nossas flags.

Uso: python make_variants.py --out <dir> [--ids id1 id2 ...]
"""
import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

CATALOG = '''"""Gerado por make_variants.py. Não editar."""
import copy
import sys

sys.path.insert(0, {harness!r})
import catalog_base as _c

TOOLCHAINS = _c.TOOLCHAINS


def carregar_ambiente():
    env = _c.carregar_ambiente()
    env["GO111MODULE"] = "off"
    _c.ENV = env
    return env


def _ours(src):
    lang, ref = next((k, v) for k, v in _c.REFERENCIAS.items() if v["src"] == src)
    return {{"lang": lang, "build": copy.deepcopy(ref["build"]), "run": list(ref["run"])}}


PROJETOS = {{
{entries}
}}

for _p in PROJETOS.values():
    _v = _p["traduzida"]
    if _v["build"] and _v["build"][0] == "cpp" and _v["run"][0].startswith("./"):
        _v["build"][1]["bin"] = _v["run"][0][2:]
'''


def read_manifest():
    with open(HERE / "variants" / "manifest.csv", newline="") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", nargs="*")
    args = ap.parse_args()

    rows = read_manifest()
    if args.ids:
        known = {r["id"] for r in rows}
        missing = [i for i in args.ids if i not in known]
        if missing:
            sys.exit(f"unknown variant(s): {', '.join(missing)}")
        rows = [r for r in rows if r["id"] in args.ids]

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    entries = []
    for r in rows:
        dest = out / r["id"]
        if dest.exists():
            shutil.rmtree(dest)
        if r["base"].startswith("ours/"):
            src = r["base"].split("/", 1)[1]
            shutil.copytree(HERE.parent / src, dest,
                            ignore=shutil.ignore_patterns("*.class", "obj", "bin", "target", ".build", "solve", "main"))
            entry = f'{{"traduzida": _ours({src!r})}}'
        else:
            shutil.copytree(HERE / "translations" / r["base"], dest)
            entry = f'{{"traduzida": copy.deepcopy(_c.PROJETOS[{r["base"]!r}]["traduzida"])}}'

        scripts = set()
        for p in r["patches"].split():
            subprocess.run(["patch", "-p1", "-s", "--no-backup-if-mismatch", "-d", str(dest),
                            "-i", str(HERE / "variants" / p)], check=True)
            scripts |= set(re.findall(r"^\+\+\+ b/(run_[^/\s]*\.sh)", (HERE / "variants" / p).read_text(), re.M))
        for script in scripts:
            with open(dest / script) as f:
                f.readline()
                cmd = [tok for tok in shlex.split(f.readline()) if tok != "$@"]
            entry = f'{{"traduzida": dict({entry[14:-1]}, run={cmd!r})}}'
        entries.append(f'    {r["id"]!r}: {entry},')

    (out / "benchmark.py").write_text(CATALOG.format(harness=str(HERE / "harness"),
                                                     entries="\n".join(entries)))
    print(f"{len(rows)} variant(s) in {out}")


if __name__ == "__main__":
    main()

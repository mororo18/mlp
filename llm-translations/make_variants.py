"""Monta as variantes do experimento e o catálogo que o bench.py -t lê.

Cada linha de variants/manifest.csv é uma variante: uma cópia de
translations/<base>/ com os patches da coluna `patches` aplicados em ordem
(caminhos relativos a variants/, separados por espaço). Nada é editado à mão.

Saída em --out: uma pasta por variante e um benchmark.py com PROJETOS,
TOOLCHAINS e carregar_ambiente(), que é a interface que o bench.py espera.
A configuração de build e execução de cada variante é a da sua base no
catálogo do harness antigo (harness/catalog_base.py), com as nossas flags.

Uso: python make_variants.py --out <dir> [--ids id1 id2 ...]
"""
import argparse
import csv
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


PROJETOS = {{
{entries}
}}
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
        shutil.copytree(HERE / "translations" / r["base"], dest)
        for p in r["patches"].split():
            subprocess.run(["patch", "-p1", "-s", "--no-backup-if-mismatch", "-d", str(dest),
                            "-i", str(HERE / "variants" / p)], check=True)
        entries.append(f'    {r["id"]!r}: {{"traduzida": copy.deepcopy(_c.PROJETOS[{r["base"]!r}]["traduzida"])}},')

    (out / "benchmark.py").write_text(CATALOG.format(harness=str(HERE / "harness"),
                                                     entries="\n".join(entries)))
    print(f"{len(rows)} variant(s) in {out}")


if __name__ == "__main__":
    main()

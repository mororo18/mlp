import re, shutil, subprocess, sys
from pathlib import Path
LT = Path('/home/victor/repos/tcc/mlp-main/llm-translations')
OURS = LT.parent / 'java'
W = Path(__import__('os').environ.get('WORK', '/tmp/variant-work')).resolve()

def mk(name, src, edit, pre=()):
    a, b = W / name / 'a', W / name / 'b'
    shutil.rmtree(W / name, ignore_errors=True)
    for d in (a, b):
        shutil.copytree(src, d, ignore=shutil.ignore_patterns('*.class', 'notas', 'CLAUDE.md', 'RESULTS.md'))
        for p in pre:
            subprocess.run(['patch', '-p1', '-s', '-d', d, '-i', LT / 'variants' / p], check=True)
    edit(b)
    out = subprocess.run(['diff', '-ruN', 'a', 'b'], cwd=W / name, capture_output=True, text=True).stdout
    assert out, name
    return out

def sub(path, pat, rep, count_min=1, flags=0):
    t = path.read_text()
    t2, n = re.subn(pat, rep, t, flags=flags)
    assert n >= count_min, (path, pat, n)
    path.write_text(t2)
    return n


def write(jobs):
    for rel, (src, f, pre) in jobs.items():
        p = LT / 'variants' / rel
        p.parent.mkdir(exist_ok=True)
        p.write_text(mk(rel.replace('/', '_'), src, f, pre))
        print('wrote', rel)

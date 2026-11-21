from common import *
OURS = LT.parent / 'java'

def z_final(b):
    sub(b / 'tSolution.java', r'private static int Z = 3;', 'private static final int Z = 3;')

def direct(b):
    sub(b / 'tSolution.java', r'private tArray s;', 'tArray s;')
    sub(b / 'tSolution.java', r'private double \[\] seq2d;', 'double [] seq2d;')
    sub(b / 'tSolution.java', r'private int n1;', 'int n1;')
    sub(b / 'tSolution.java', r'private static int Z = 3;', 'static int Z = 3;')
    sub(b / 'tArray.java', r'private int \[\] arr;', 'int [] arr;')
    sub(b / 'tInfo.java', r'private double  \[\]\[\] c;', 'double  [][] c;')
    g = b / 'GILS_RVND.java'
    t = g.read_text()
    n0 = t.count('getSeq(') + t.count('setSeq(') + t.count('getCost(') + t.count('getSolut().get(')
    t = re.sub(r'solut\.getSolut\(\)\.get\(([^()]+)\)', r'solut.s.arr[\1]', t)
    t = re.sub(r'\bs\.get\(([^()]+)\)', r's.arr[\1]', t)
    t = re.sub(r'solut\.setSeq\(([^,()]+),([^,()]+),\s*info\.([TCW]),\s*(.*)\);',
               r'solut.seq2d[((\1)*solut.n1 + (\2))*tSolution.Z + tInfo.\3] = \4;', t)
    t = re.sub(r'solut\.getSeq\(([^,()]+),([^,()]+),\s*info\.([TCW])\)',
               r'solut.seq2d[((\1)*solut.n1 + (\2))*tSolution.Z + tInfo.\3]', t)
    t = re.sub(r'info\.getCost\(([^()]+?),\s*([^()]+?)\)', r'info.c[\1][\2]', t)
    left = re.findall(r'setSeq\(|info\.getCost\(|getSolut\(\)\.get\(', t)
    assert not left and t.count('getSeq(') == 1, left
    g.write_text(t)
    print('direct: replaced', n0)

def no_locals(b):
    g = b / 'GILS_RVND.java'
    out, m, n = [], {}, 0
    for line in g.read_text().splitlines(keepends=True):
        d = re.match(r'\s*int (s_\w+) = (solut\.getSolut\(\)\.get\([^;]*\));\s*$', line)
        if d:
            m[d.group(1)] = d.group(2); n += 1
            continue
        for k in sorted(m, key=len, reverse=True):
            line = re.sub(r'\b%s\b' % k, m[k], line)
        out.append(line)
    g.write_text(''.join(out))
    assert not re.search(r'\bs_\w+\b', g.read_text())
    print('no_locals: removed', n)

def arraycopy(b):
    sub(b / 'tArray.java', r'    private void shift\(int from, int to, int sz\) \{.*?\n    \}\n',
        '    private void shift(int from, int to, int sz) {\n        System.arraycopy(arr, from, arr, to, sz);\n    }\n',
        flags=re.S)

def aos(b):
    m = b / 'Main.java'
    sub(m, r'        double\[\] seqT;\n        double\[\] seqC;\n        double\[\] seqW;\n', '        double[] seq;\n')
    sub(m, r'        solut\.seqT = new double\[sz\];\n        solut\.seqC = new double\[sz\];\n        solut\.seqW = new double\[sz\];\n',
        '        solut.seq = new double[3 * sz];\n')
    for f, o in (('T', 0), ('C', 1), ('W', 2)):
        sub(m, r'\.seq%s\[to_1D\(([^()]+)\)\]' % f, r'.seq[3 * to_1D(\1) + %d]' % o)
    assert not re.search(r'seq[TCW]', m.read_text())

jobs = {
    'ours_java/z_final.patch': (OURS, z_final, ()),
    'ours_java/direct_access.patch': (OURS, direct, ()),
    'ours_java/no_local_copies.patch': (OURS, no_locals, ()),
    'ours_java/arraycopy_shift.patch': (OURS, arraycopy, ()),
    'cpp_java/aos.patch': (LT / 'translations/cpp_java', aos, ('cpp_java/measure.patch',)),
}
write(jobs)

from common import *
OURS = LT.parent / 'swift'

def direct(b):
    m = b / 'Sources/ILS.swift'
    n = sub(m, r's\.getSeq\(([^()]+),\s*([^()]+)\)', r's.seq[(\1) * s.dimension + (\2)]')
    k = sub(m, r'data\.getDistance\(([^()]+?),\s*([^()]+?)\)', r'data.matrix[\1][\2]')
    t = m.read_text(); assert 'getSeq(' not in t and 'getDistance(' not in t
    m2 = b / 'Sources/Solution.swift'
    a = sub(m2, r'self\.getSeq\(([^()]+?),\s*([^()]+?)\)', r'seq[(\1) * dimension + (\2)]')
    c = sub(m2, r'data\.getDistance\(([^()]+?),\s*([^()]+?)\)', r'data.matrix[\1][\2]')
    print('direct', n, k, a, c)

def no_locals(b):
    m = b / 'Sources/ILS.swift'
    out, mp, n = [], {}, 0
    for line in m.read_text().splitlines(keepends=True):
        d = re.match(r'\s*let (s[A-Z]\w*) = (s\.sequence\[[^\]]*\])\s*$', line)
        if d:
            mp[d.group(1)] = d.group(2); n += 1
            continue
        for k in sorted(mp, key=len, reverse=True):
            line = re.sub(r'\b%s\b' % k, mp[k], line)
        out.append(line)
    m.write_text(''.join(out))
    print('no_locals', n)

def route_copy(b):
    m = b / 'Sources/ILS.swift'
    sub(m, r'func perturb\(s: Solution, data: inout Data\) -> Solution \{\n\tvar s = s\n',
        'func perturb(s src: Solution, into s: inout Solution, data: inout Data) {\n\ts.sequence = src.sequence\n')
    sub(m, r'\ts\.updateSubseqMatrix\(data: data\)\n\n\treturn s\n\}\n\nfunc searchSwap',
        '\ts.updateSubseqMatrix(data: data)\n}\n\nfunc searchSwap')
    sub(m, r'\t\t\t\tlocalBest = current\n', '\t\t\t\tlocalBest.sequence = current.sequence\n\t\t\t\tlocalBest.cost = current.cost\n')
    sub(m, r'current = perturb\(s: localBest, data: &data\)', 'perturb(s: localBest, into: &current, data: &data)')

write({'ours_swift/direct_access.patch': (OURS, direct, ()),
       'ours_swift/no_local_copies.patch': (OURS, no_locals, ()),
       'ours_swift/route_copy.patch': (OURS, route_copy, ())})

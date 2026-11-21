from common import *
OURS = LT.parent / 'go'

def eps(b):
    m = b / 'main.go'
    sub(m, r'type tSubseq struct', 'const dblEpsilon = 2.220446049250313e-16\n\ntype tSubseq struct')
    n = sub(m, r'cost_best = cost_new\n', 'cost_best = cost_new - dblEpsilon\n')
    sub(m, r'if cost_best < solut\.at\(0, data\.dimen\)\.C \{', 'if cost_best < solut.at(0, data.dimen).C-dblEpsilon {')
    k = sub(m, r'if cost_best < solut\.cost \{', 'if cost_best < solut.cost-dblEpsilon {')
    print('eps', n, k)

def direct(b):
    m = b / 'main.go'
    n = sub(m, r'solut\.at\(([^()]+),\s*([^()]+)\)', r'solut.seq[(\1)*solut.dim+(\2)]')
    assert '.at(' not in m.read_text()
    print('direct', n)

write({'ours_go/epsilon.patch': (OURS, eps, ()),
       'ours_go/direct_index.patch': (OURS, direct, ())})

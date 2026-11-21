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

def style(b):
    m = b / 'main.go'
    sub(m, r'\tT float64\n\tW float64\n\tC float64\n', '\tT float64\n\tC float64\n\tW float64\n')
    sub(m, r'\tn_list_b := \[\]int\{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3\}\n\tn_list := make\(\[\]int, 5\)\n\n\tcopy\(n_list, n_list_b\)\n',
        '\tn_list := []int{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}\n')
    sub(m, r'\t\t\tn_list = make\(\[\]int, 5\)\n\t\t\tcopy\(n_list, n_list_b\)\n',
        '\t\t\tn_list = []int{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}\n')
    sub(m, r'n_list = remove\(n_list, index\)', 'n_list = append(n_list[:index], n_list[index+1:]...)')

write({'ours_go/style.patch': (OURS, style, ('ours_go/direct_index.patch',)),'ours_go/epsilon.patch': (OURS, eps, ()),
       'ours_go/direct_index.patch': (OURS, direct, ())})

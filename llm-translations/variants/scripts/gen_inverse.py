from common import *
T = LT / 'translations'

def go_at(b):
    m = b / 'go_src' / 'main.go'
    sub(m, r'type tSolution struct \{\n\ts    \[\]int\n\tseq  \[\]tSeqInfo\n\tcost float64\n\}\n',
        'type tSolution struct {\n\ts    []int\n\tseq  []tSeqInfo\n\tdim  int\n\tcost float64\n}\n\n'
        'func (s *tSolution) at(i int, j int) *tSeqInfo {\n\treturn &s.seq[i*s.dim+j]\n}\n')
    sub(m, r'(\tsolut\.cost = dblMax\n\n\treturn solut\n)', r'\tsolut.dim = data.dimen + 1\n\1')
    n = sub(m, r'solut\.seq\[to_1D\(([^()]+), ([^()]+), data\.dimen\)\]', r'solut.at(\1, \2)')
    assert 'solut.seq[to_1D' not in m.read_text()
    print('go_at', n)

def swift_full(b):
    m = b / 'main.swift'
    sub(m, r'func solutionCpy\(_ src: TSolution, _ tgt: inout TSolution\) \{\n    tgt\.s = src\.s\n    tgt\.cost = src\.cost\n\}\n',
        'func solutionCpy(_ src: TSolution, _ tgt: inout TSolution) {\n    tgt = src\n}\n')
    sub(m, r'(\n            )solutCrnt\.s = perturb\(solutPartial, &data\)',
        r'\1solutCrnt = solutPartial\1solutCrnt.s = perturb(solutPartial, &data)')

write({'cpp_go/at_accessor.patch': (T / 'cpp_go', go_at, ('cpp_go/measure.patch',)),
       'cpp_swift/full_copy.patch': (T / 'cpp_swift', swift_full, ('cpp_swift/measure.patch',))})

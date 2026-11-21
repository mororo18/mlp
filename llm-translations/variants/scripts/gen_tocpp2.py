from common import *
T = LT / 'translations'

def sw_route(b):
    m = b / 'src/ILS.cpp'
    sub(m, r'Solution perturb\(Solution s, Data &data\) \{\n', 'void perturb(const Solution &src, Solution &s, Data &data) {\n\ts.sequence = src.sequence;\n')
    t = m.read_text(); i = t.index('void perturb(')
    j = t.index('\ts.updateSubseqMatrix(data);\n\n\treturn s;\n}\n', i)
    t = t[:j] + '\ts.updateSubseqMatrix(data);\n}\n' + t[j + len('\ts.updateSubseqMatrix(data);\n\n\treturn s;\n}\n'):]
    m.write_text(t)
    sub(m, r'([ \t]+)localBest = current;\n', r'\1localBest.sequence = current.sequence;\n\1localBest.cost = current.cost;\n')
    sub(m, r'current = perturb\(localBest, data\);', 'perturb(localBest, current, data);')

def sw_idx(b):
    sub(b / 'src/Solution.hpp', r'return seq\[i \* dimension \+ j\];', 'return seq[(unsigned)(i * dimension + j)];')
    sub(b / 'src/Data.hpp', r'return matrix\[i \* dimension \+ j\];', 'return matrix[(unsigned)(i * dimension + j)];')

def sw_rows(b):
    sub(b / 'src/Data.hpp', r'std::vector<double> matrix;', 'double **matrix = nullptr;')
    sub(b / 'src/Data.hpp', r'return matrix\[i \* dimension \+ j\];', 'return matrix[i][j];')
    d = b / 'src/Data.cpp'
    sub(d, r'this->matrix = std::vector<double>\(dimension \* dimension, 0\.0\);',
        'this->matrix = new double *[dimension];\n\tfor (int r = 0; r < dimension; r++)\n\t\tthis->matrix[r] = new double[dimension]();')
    n = sub(d, r'this->matrix\[(\w+) \* dimension \+ (\w+)\]', r'this->matrix[\1][\2]')
    assert 'dimension +' not in re.sub(r'.*\n', lambda x: x.group(0) if 'matrix' in x.group(0) else '', d.read_text())
    print('sw rows', n)

def cs_idx(b):
    n = sub(b / 'tSolution.hpp', r'seq\[i \* \(d \+ 1\) \* 3 \+ j \* 3 \+ k\]', 'seq[(unsigned)(i * (d + 1) * 3 + j * 3 + k)]')
    sub(b / 'Data.hpp', r'return c\[to_1D\(i, j\)\];', 'return c[(unsigned)to_1D(i, j)];')
    print('cs idx', n)

def cs_rows(b):
    h = b / 'Data.hpp'
    sub(h, r'        std::vector<double> c;\n', '        double **c = nullptr;\n')
    sub(h, r'return c\[to_1D\(i, j\)\];', 'return c[i][j];')
    d = b / 'Data.cpp'
    sub(d, r'        c\.assign\(\(size_t\)dimension \* dimension, 0\.0\);\n',
        '        c = new double *[dimension];\n        for (int r = 0; r < dimension; r++)\n            c[r] = new double[dimension]();\n')
    sub(d, r'c\[to_1D\((\w), (\w)\)\] = ', r'c[\1][\2] = ', count_min=2)

write({'swift_cpp/route_copy.patch': (T / 'swift_cpp', sw_route, ('swift_cpp/measure.patch',)),
       'swift_cpp/index32.patch': (T / 'swift_cpp', sw_idx, ('swift_cpp/measure.patch',)),
       'swift_cpp/row_matrix.patch': (T / 'swift_cpp', sw_rows, ('swift_cpp/measure.patch',)),
       'csharp_cpp/index32.patch': (T / 'csharp_cpp', cs_idx, ('csharp_cpp/measure.patch',)),
       'csharp_cpp/row_matrix.patch': (T / 'csharp_cpp', cs_rows, ('csharp_cpp/measure.patch',))})

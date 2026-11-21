from common import *
T = LT / 'translations'
EPS = 'std::numeric_limits<double>::epsilon()'

def go_eps(b):
    m = b / 'main.cpp'
    n = sub(m, r'cost_best = cost_new;', 'cost_best = cost_new - %s;' % EPS)
    sub(m, r'if \(cost_best < solut->at\(0, data->dimen\)->C\) \{', 'if (cost_best < solut->at(0, data->dimen)->C - %s) {' % EPS)
    k = sub(m, r'if \(cost_best < solut->cost\) \{', 'if (cost_best < solut->cost - %s) {' % EPS)
    t = m.read_text()
    if '#include <limits>' not in t:
        m.write_text('#include <limits>\n' + t)
    print('go eps', n, k)

def go_idx(b):
    sub(b / 'main.cpp', r'tSubseq \*at\(int i, int j\) \{ return &seq\[i \* dim \+ j\]; \}',
        'tSubseq *at(int i, int j) { return &seq[(unsigned)(i * dim + j)]; }')

def go_rows(b):
    m = b / 'main.cpp'
    sub(m, r'struct tData \{\n    std::vector<std::vector<double>> c;', 'struct tData {\n    double **c;')
    sub(m, r'    read_data\(verbose, data\.dimen, data\.c, data\.rnd\);\n',
        '    std::vector<std::vector<double>> cost_vec;\n'
        '    read_data(verbose, data.dimen, cost_vec, data.rnd);\n'
        '    data.c = new double *[data.dimen];\n'
        '    for (int i = 0; i < data.dimen; i++) {\n'
        '        data.c[i] = new double[data.dimen];\n'
        '        for (int j = 0; j < data.dimen; j++)\n'
        '            data.c[i][j] = cost_vec[i][j];\n'
        '    }\n')

def f_eps(b):
    m = b / 'main.cpp'
    n = sub(m, r'cost_best = cost_new;', 'cost_best = cost_new - %s;' % EPS)
    k = sub(m, r'if \(cost_best < solut\.cost\) \{', 'if (cost_best < solut.cost - %s) {' % EPS)
    t = m.read_text()
    if '#include <limits>' not in t:
        m.write_text('#include <limits>\n' + t)
    print('fortran eps', n, k)

def f_idx(b):
    n = sub(b / 'types.hpp',
            r'return seq\[static_cast<std::size_t>\(i - 1\) \* static_cast<std::size_t>\(seq_dim\) \+ static_cast<std::size_t>\(j - 1\)\];',
            'return seq[(unsigned)((i - 1) * seq_dim + (j - 1))];')
    print('fortran idx', n)

def f_rows(b):
    h = b / 'types.hpp'
    sub(h, r'struct tMatrix \{\n    int n = 0;\n    std::vector<typeReal> data;\n\n    void allocate\(int dimen\) \{\n        n = dimen;\n        data\.assign\([^\n]*\n    \}\n',
        'struct tMatrix {\n    int n = 0;\n    typeReal **data = nullptr;\n\n    void allocate(int dimen) {\n        n = dimen;\n'
        '        data = new typeReal *[dimen];\n        for (int i = 0; i < dimen; i++)\n            data[i] = new typeReal[dimen]();\n    }\n')
    k = sub(h, r'return data\[static_cast<std::size_t>\(j - 1\) \* static_cast<std::size_t>\(n\) \+ static_cast<std::size_t>\(i - 1\)\];',
            'return data[i - 1][j - 1];')
    print('fortran rows', k)

write({'go_cpp/epsilon.patch': (T / 'go_cpp', go_eps, ('go_cpp/measure.patch',)),
       'go_cpp/index32.patch': (T / 'go_cpp', go_idx, ('go_cpp/measure.patch',)),
       'go_cpp/row_matrix.patch': (T / 'go_cpp', go_rows, ('go_cpp/measure.patch',)),
       'fortran_cpp/epsilon.patch': (T / 'fortran_cpp', f_eps, ()),
       'fortran_cpp/index32.patch': (T / 'fortran_cpp', f_idx, ()),
       'fortran_cpp/row_matrix.patch': (T / 'fortran_cpp', f_rows, ())})

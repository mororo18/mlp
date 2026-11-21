from common import *
OURS = LT.parent / 'cplusplus'
TR = LT / 'translations/rust_cpp'

def idx64(b):
    sub(b / 'main.cpp', r'inline unsigned to_1D\(int x, int y, int size\) \{\n    int a = x\*\(size\+1\) \+ y;\n',
        'inline std::size_t to_1D(std::size_t x, std::size_t y, std::size_t size) {\n    std::size_t a = x*(size+1) + y;\n')

def idx32(b):
    h = b / 'cpp/subseq.hpp'
    n = sub(h, r'return data\[row_size \* row \+ col\];', 'return data[(unsigned)(row_size * row + col)];')
    print('idx32', n)

def rows(b):
    h = b / 'cpp/subseq.hpp'
    t = h.read_text()
    i = t.index('class CostMatrix')
    c = t[i:]
    c = c.replace('#ifdef FLAT\n    std::size_t row_size;\n    std::vector<double> data;\n', '#ifdef FLAT\n    double **data = nullptr;\n', 1)
    c = c.replace('#ifdef FLAT\n        : row_size(n), data(n * n, 0.0) {\n    }\n',
                  '#ifdef FLAT\n    {\n        data = new double *[n];\n        for (std::size_t i = 0; i < n; i++)\n            data[i] = new double[n]();\n    }\n', 1)
    k = c.count('return data[row_size * row + col];')
    c = c.replace('return data[row_size * row + col];', 'return data[row][col];')
    assert k == 2, k
    h.write_text(t[:i] + c)

def route_int(b):
    h = b / 'cpp/subseq.hpp'
    sub(h, r'    std::vector<std::size_t> s;\n', '    std::vector<int> s;\n')

write({'ours_cplusplus/index64.patch': (OURS, idx64, ()),
       'rust_cpp/index32.patch': (TR, idx32, ('rust_cpp/measure.patch',)),
       'rust_cpp/row_matrix.patch': (TR, rows, ('rust_cpp/measure.patch',)),
       'rust_cpp/route_int.patch': (TR, route_int, ('rust_cpp/measure.patch',))})

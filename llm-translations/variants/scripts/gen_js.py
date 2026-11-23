from common import *
OURS = LT.parent / 'javascript'

def typed(b):
    sub(b / 'main.js', r'    seq = Array\(\(dimension\+1\)\*\(dimension\+1\)\*3\);', '    seq = new Float64Array((dimension+1)*(dimension+1)*3);')

def consts(b):
    m = b / 'main.js'
    t = m.read_text()
    n = 0
    for f, v in (('T', 0), ('C', 1), ('W', 2)):
        t, k = re.subn(r'info\.%s\b' % f, str(v), t); n += k
    m.write_text(t); print('consts', n)

def locals_(b):
    sub(b / 'main.js', r'    seq = Array\(', '    var seq = Array(')
    sub(b / 'main.js', r'    subseq = subseq_fill\(', '    var subseq = subseq_fill(')

def fdecl(b):
    sub(b / 'main.js', r'var to_1D = function \(x, y, z, n\) \{\n    return \(x\*n \+ y\)\*3 \+ z;\n\}\n',
        'function to_1D(x, y, z, n) {\n    return (x*n + y)*3 + z;\n}\n')

def tr_array(b):
    for f in 'TCW':
        sub(b / 'main.js', r'    solut\.seq_%s = new Float64Array\(sz\);' % f, '    solut.seq_%s = Array(sz);' % f)

def three_arrays(b, mk='new Float64Array(sz)'):
    # P20: T, C e W em três Float64Array, como na tradução C++ -> JavaScript
    m = b / 'main.js'
    sub(m, r'    seq = Array\(\(dimension\+1\)\*\(dimension\+1\)\*3\);',
        '    var sz = (dimension+1)*(dimension+1);\n    seq = [%s, %s, %s];' % (mk, mk, mk))
    sub(m, r'var to_1D = function \(x, y, z, n\) \{\n    return \(x\*n \+ y\)\*3 \+ z;\n\}\n',
        'var to_2D = function (x, y, n) {\n    return x*n + y;\n}\n')
    sub(m, r'\bsubseq\[to_1D\(([^,()]+), ([^,()]+), info\.([TCW]), info\.dimension\)\]', r'subseq[info.\3][to_2D(\1, \2, info.dimension)]', 3)
    n = sub(m, r'\bseq\[to_1D\(([^,()]+), ([^,()]+), info\.([TCW]), info\.dimension\)\]', r'seq_\3[to_2D(\1, \2, info.dimension)]', 50)
    k = sub(m, r'(\nfunction (subseq_load|search_swap|search_two_opt|search_reinsertion)\(s, seq, info(, opt)?\) \{\n)',
            r'\1    var seq_T = seq[info.T], seq_C = seq[info.C], seq_W = seq[info.W];\n', 4)
    assert k == 4 and 'to_1D' not in m.read_text(), (n, k)
    print('three_arrays', n)

def tr_interleaved(b, mk='new Float64Array(sz*3)'):
    # P20: uma Float64Array intercalada (T, C, W), como na nossa implementação
    m = b / 'main.js'
    sub(m, r'    solut\.seq_T = new Float64Array\(sz\);\n    solut\.seq_C = new Float64Array\(sz\);\n    solut\.seq_W = new Float64Array\(sz\);',
        '    solut.seq = %s;' % mk)
    k = sub(m, r'const T = solut\.seq_T, C = solut\.seq_C, W = solut\.seq_W;', 'const S = solut.seq;', 4)
    off = {'T': '', 'C': ' + 1', 'W': ' + 2'}
    n = sub(m, r'\b([TCW])\[to_1D\(([^()]*)\)\]', lambda mo: 'S[to_1D(%s)*3%s]' % (mo.group(2), off[mo.group(1)]), 30)
    assert 'solut.seq_' not in m.read_text() and k == 4, k
    print('tr_interleaved', n)

def three_arrays_array(b):
    # P20 sobre o Array: T, C e W em três Array, como a tradução com o tr_array
    three_arrays(b, 'Array(sz)')

def tr_interleaved_array(b):
    # P20 sobre o Array: um Array intercalado, como na nossa implementação
    tr_interleaved(b, 'Array(sz*3)')

HOT = ('subseq_load', 'search_swap', 'search_two_opt', 'search_reinsertion')

def local_info(b):
    # P21: info.dimension e info.c em variáveis locais no início das buscas, como
    # a tradução C++ -> JavaScript faz com dimen e cost
    m = b / 'main.js'
    parts = re.split(r'(?m)^(?=function )', m.read_text())
    k = 0
    for i, p in enumerate(parts):
        mo = re.match(r'function (\w+)\([^)]*\) \{\n', p)
        if not mo or mo.group(1) not in HOT:
            continue
        body = p[mo.end():]
        body = re.sub(r'info\.dimension\b', 'dimension', body)
        body = re.sub(r'info\.c\[', 'c[', body)
        parts[i] = mo.group(0) + '    var dimension = info.dimension, c = info.c;\n' + body
        k += 1
    assert k == 4, k
    m.write_text(''.join(parts))

def local_info_consts(b):
    local_info(b)
    consts(b)

write({'cpp_javascript/array.patch': (LT / 'translations' / 'cpp_javascript', tr_array, ('cpp_javascript/measure.patch',)),
       'cpp_javascript/interleaved.patch': (LT / 'translations' / 'cpp_javascript', tr_interleaved, ('cpp_javascript/measure.patch',)),
       'cpp_javascript/interleaved_array.patch': (LT / 'translations' / 'cpp_javascript', tr_interleaved_array, ('cpp_javascript/measure.patch',)),
       'ours_javascript/three_arrays_array.patch': (OURS, three_arrays_array, ()),
       'ours_javascript/float64array.patch': (OURS, typed, ()),
       'ours_javascript/three_arrays.patch': (OURS, three_arrays, ()),
       'ours_javascript/const_offsets.patch': (OURS, consts, ()),
       'ours_javascript/no_globals.patch': (OURS, locals_, ()),
       'ours_javascript/function_decl.patch': (OURS, fdecl, ()),
       'ours_javascript/local_info.patch': (OURS, local_info, ()),
       'ours_javascript/local_info_consts.patch': (OURS, local_info_consts, ())})

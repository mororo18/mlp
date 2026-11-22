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

def ncr(b):
    sub(b / 'run_node.sh', r'node main\.js', 'node --no-concurrent-recompilation main.js')

write({'ours_javascript/no_concurrent_recompilation.patch': (OURS, ncr, ()),
       'cpp_javascript/array.patch': (LT / 'translations' / 'cpp_javascript', tr_array, ('cpp_javascript/measure.patch',)),'ours_javascript/float64array.patch': (OURS, typed, ()),
       'ours_javascript/const_offsets.patch': (OURS, consts, ()),
       'ours_javascript/no_globals.patch': (OURS, locals_, ()),
       'ours_javascript/function_decl.patch': (OURS, fdecl, ())})

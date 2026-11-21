from common import *
OURS = LT.parent / 'julia'

def blocks(t, names, fn):
    for name in names:
        m = re.search(r'\nfunction %s\(.*?\n(.*?)\nend\n' % name, t, re.S)
        assert m, name
        t = t[:m.start(1)] + fn(m.group(1)) + t[m.end(1):]
    return t

def hoist(b):
    m = b / 'main.jl'
    t = m.read_text()
    t = t.replace('const C = 3\n', 'const C = 3\nconst DBL_EPSILON = eps(Float64)\n', 1)
    n = t.count('data.EPSILON'); t = t.replace('data.EPSILON', 'DBL_EPSILON')
    def f(body):
        return '    dimen = data.dimen\n' + body.replace('data.dimen', 'dimen')
    t = blocks(t, ['update_subseq_info_matrix', 'search_swap', 'search_two_opt', 'search_reinsertion'], f)
    m.write_text(t); print('hoist eps', n)

def structs(b):
    m = b / 'main.jl'
    t = m.read_text()
    t = t.replace('mutable struct tSolution\n', 'struct tSI\n    T::Float64\n    C::Float64\n    W::Float64\nend\n\nmutable struct tSolution\n', 1)
    t = t.replace('    seq::Array{Float64, 3}\n', '    seq::Vector{tSI}\n', 1)
    k = t.count('Array{Float64, 3}(undef, 3, data.dimen+1, data.dimen+1)')
    t = t.replace('Array{Float64, 3}(undef, 3, data.dimen+1, data.dimen+1)', 'Vector{tSI}(undef, (data.dimen+1)*(data.dimen+1))')
    upd = '''    n = data.dimen+1
    for i in 1:n
        k::Int = 1 - i - (i==1)
        base = (i-1)*n
        solut.seq[base + i] = tSI(0.0, 0.0, i!=1)
        for j in i+1:n
            j_prev = j-1
            prev = solut.seq[base + j_prev]
            tt = data.cost[solut.s[j_prev], solut.s[j]] + prev.T
            solut.seq[base + j] = tSI(tt, tt + prev.C, j + k)
        end
    end
    solut.cost = solut.seq[n].C'''
    t = blocks(t, ['update_subseq_info_matrix'], lambda body: upd)
    t, n = re.subn(r'solut\.seq\[([TCW]), ([^,\]]+), ([^\]]+)\]', r'solut.seq[((\3)-1)*(data.dimen+1) + (\2)].\1', t)
    assert 'solut.seq[T' not in t and 'solut.seq[C' not in t and 'solut.seq[W' not in t
    m.write_text(t); print('structs', k, n)

def zero_based(b):
    m = b / 'main.jl'
    sub(m, r'\n    push!\(s, 1\)\n    return s\n', '\n    push!(s, 1)\n    return s .- 1\n')
    n = sub(m, r'data\.cost\[solut\.s\[([^\[\]]+)\], solut\.s\[([^\[\]]+)\]\]', r'data.cost[solut.s[\1]+1, solut.s[\2]+1]')
    assert not re.search(r'data\.cost\[solut\.s\[[^\]]+\], ', m.read_text())
    print('zero_based', n)

write({'ours_julia/zero_based.patch': (OURS, zero_based, ()),
       'ours_julia/hoist_fields.patch': (OURS, hoist, ()),
       'ours_julia/struct_vector.patch': (OURS, structs, ())})

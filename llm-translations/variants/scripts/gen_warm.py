from common import *

def ours_warm(b):
    sub(b / 'main.java', r'(    private static void read\(boolean verbose\)\{\n)',
        r'\1        new GILS_RVND(verbose).solve();\n')

def cpp_warm(b):
    sub(b / 'Main.java', r'(\n        long t1 = System\.nanoTime\(\);)',
        r'\n        GILS_RVND(Imax, Iils, data, verbose);\n        data.rnd_index = 0;\1')

def ours_julia_warm(b):
    sub(b / 'main.jl', r'(\n    _cpu0 = _cputime\(\))',
        r'\n    GILS_RVND(Imax, Iils, R, data, verbose)\n    data.rnd_index = 1\1')

def cpp_julia_warm(b):
    sub(b / 'main.jl', r'(\n    t1 = time_ns\(\))',
        r'\n    GILS_RVND(Imax, Iils, data, verbose)\n    data.rnd_index = 1\1')

write({
    'ours_julia/warm.patch': (LT.parent / 'julia', ours_julia_warm, ()),
    'cpp_julia/warm.patch': (LT / 'translations' / 'cpp_julia', cpp_julia_warm, ('cpp_julia/measure.patch',)),
    'ours_java/warm.patch': (OURS, ours_warm, ()),
    'cpp_java/warm.patch': (LT / 'translations' / 'cpp_java', cpp_warm, ('cpp_java/measure.patch',)),
})

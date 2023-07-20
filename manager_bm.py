import os
import argparse
import time

bm_dir = '../virtu_test'
inst_list_file = 'virtu_test.txt' 
inst_list = []
min_test = 1

#sources = ["java", "dotnet", "python3", "pypy", "julia", "cpp", "cpp-OOP",
         #"fortran", "node", "luajit", "rust", "c", "golang"]
sources = ["java", "dotnet", "python3", "pypy", "julia", "cpp", "cpp-OOP",
         "fortran", "node", "lua", "luajit", "rust", "c", "matlab", "golang"]

lang_dir = {
        "dotnet": "csharp",
        "mcs": "csharp",
        "java": "java",
        "python3": "python",
        "pypy": "python",
        "julia": "julia",
        "cpp" : "cplusplus",
        "cpp-OOP" : "cppOOP",
        "c" : "c",
        "fortran" : "fortran",
        "node" : "javascript",
        "lua" : "lua",
        "luajit" : "lua",
        "rust" : "rust",
        "matlab" : "octave",
        "golang" : "go"
        }


if os.path.isdir(bm_dir) == False:
    os.mkdir(bm_dir)

with open(inst_list_file) as f:
    for line in f:
        inst_list.append(line.replace('\n', ''))

def count(source, lang, inst, path):

    print(bm_dir)
    f_name = os.path.join(bm_dir, lang+'.csv')
    c = 0
    if os.path.isfile(f_name) == False:
        return c
    with open(f_name) as f:
        for line in f:
            if line.find(source) >= 0 and line.find(inst) >= 0:
                c += 1

    return c


parser = argparse.ArgumentParser(description='Rodador Tudao')
parser.add_argument('--lang' , nargs='+', required=True, help='Sources: python3, java, mcs, dotnet, julia, cpp, lua, javascript, matlab, golang')
parser.add_argument('--min' ,  default=min_test, type=int, help='Quantidade minima de rodadas de cada linguagem')
parser.add_argument('--out' , default=bm_dir, help='Output dir')

args = parser.parse_args()

min_test = args.min
bm_dir = args.out

for i in args.lang:
    if i not in sources and i != 'all':
        print("{} is not suported".format(i))
        exit(0)

if args.lang[0] != 'all':
    sources = args.lang[:]

print(sources)
print(inst_list)
for i in range(min_test):
    for s in sources:
        for inst in inst_list:

            crnt_time = time.localtime()
            hour = crnt_time.tm_hour
            day = crnt_time.tm_wday
            _min = crnt_time.tm_min

            while (_min < 3 or _min > 5):# and day < 5:
            #while (hour >= 9 and hour < 21):# and day < 5:
                time.sleep(5)
                hour = time.localtime().tm_hour
                day = crnt_time.tm_wday
                _min = crnt_time.tm_min

            if count(s, lang_dir[s], inst, bm_dir) < min_test:
                os.system(f'python3.8 run_bm.py -i {inst} --lang {s}')





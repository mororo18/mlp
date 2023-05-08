import os
import time

bm_dir = '../mlp_test'
inst_list_file = 'inst-list'
inst_list = []

min_test = 5

#sources = ["java", "dotnet", "mcs", "python3", "pypy", "julia", "cpp", "cpp-OOP",
        #"fortran", "node", "lua", "luajit", "rust", "octave", "c", "matlab", "golang"]
sources = ["julia", "cpp",
        "fortran","c"]

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

    f_name = os.path.join(bm_dir, lang+'.csv')
    c = 0
    if os.path.isfile(f_name) == False:
        return c
    with open(f_name) as f:
        for line in f:
            if line.find(source) > 0 and line.find(inst) > 0:
                c += 1

    return c

for s in sources:
    for inst in inst_list:

        crnt_time = time.localtime()
        hour = crnt_time.tm_hour

        while hour >= 9 and hour < 21:
            time.sleep(5)
            hour = time.localtime().tm_hour

        if count(s, lang_dir[s], inst, bm_dir) < min_test:
            os.system(f'python3.8 run_bm.py -i {inst} --lang {s}')
        else:
            continue

        


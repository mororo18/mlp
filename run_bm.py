import subprocess
import argparse
import tempfile
import os
import sys
import pandas as pd
import concurrent.futures
import psutil as pU
import time

data_dir = "../mlp_testao"
###data_dir = "../mlp_tudao"

DEFAULT_TIMEOUT = 900  # segundos por execucao individual; ver doc/tcc/EXPERIMENTO.md secao 4

def get_branch():
    out = subprocess.check_output(['git', 'branch']).decode('utf-8').split('\n')
    for b in out:
        if b.find('*') != -1:
            branch = b.replace('*', '').replace(' ', '')
            return branch

def ds_open(lang_name, data_dir):
    ds_path = os.path.join(data_dir, lang_name + ".csv")

    if os.path.isfile(ds_path):
        df = pd.read_csv(ds_path)
    else:
        os.system('cp data/template.csv ' + ds_path)
        df = pd.read_csv(ds_path)

    return df

def get_mem_avg(pid, timeout):

    p = pU.Process(pid=pid)
    start = time.time()
    MB = 1024**2

    gap = 1e-10

    mem_max = -99999
    mem_avg = 0
    count = 0
    timed_out = False
    while True:
        if time.time() - start >= timeout:
            timed_out = True
            print("[timeout] estourou {}s, matando processo".format(timeout))
            try:
                p.kill()
            except pU.NoSuchProcess:
                pass
            break

        with p.oneshot():
            #if p.status() != pU.STATUS_RUNNING:
            if p.status() == pU.STATUS_ZOMBIE:
            # if p.status() == pU.STATUS_STOPPED:
                print(p.status())
                break
            #if p.status() == pU.STATUS_SLEEPING:
            #if p.status() == pU.STATUS_RUNNING:
            #if p.status() == pU.STATUS_ZOMBIE:
                #print(p.status())

            mem = p.memory_info()
            mem_total = (float(mem.rss)- float(mem.shared))/ MB
            if mem_total < gap:
                print("break by gap")
                break;
            elif mem_total > mem_max:
                mem_max = mem_total

            mem_avg += mem_total
        count += 1
        time.sleep(0.05)

    mem_avg = mem_avg / count if count else float('nan')
    mem_max = mem_max if count else float('nan')

    return {"mem_avg" : mem_avg, "mem_max" : mem_max, "mem_lookups" : count, "timed_out" : timed_out}

def get_COST(_str):
    lines = _str.split("\n")

    for l in lines :
        #print(l)
        if l.find("COST") != -1:
            i = l.find(":")
            col = l[:i]
            val = float(l[i+1:].replace(',', '.'))
            return val

def get_TIME(_str):
    lines = _str.split("\n")

    for l in lines :
        #print(l)
        if l.find("TIME") != -1:
            i = l.find(":")
            col = l[:i]
            val = float(l[i+1:].replace(',', '.'))
            return val
    return

def get_info(lang, timeout):
    f = open('run_' + lang + '.sh', 'r')
    #print(f.readline())
    f.readline()
    cmd_line = f.readline().replace('\n', '').split(' ')
    print(cmd_line)
    f.close()

    with tempfile.TemporaryFile() as tempf:
        proc = subprocess.Popen(cmd_line, stdout=tempf)
        pid = proc.pid

        info = get_mem_avg(pid, timeout)

        proc.wait()
        tempf.seek(0)
        output = str(tempf.read().decode("utf-8"))

        if info["timed_out"]:
            # saida pode estar incompleta/vazia - nao confiar em COST/TIME parciais
            COST = None
            TIME = None
        else:
            COST = get_COST(output)
            TIME = get_TIME(output)
        info.update({"COST" : [COST], "TIME" : [TIME]})
        print(info)

        return info

def main():
    global data_dir

    parser = argparse.ArgumentParser(description='Run Benchmark')
    parser.add_argument('-i' ,'--instance', help='Path to the instance file.', required= not ('--instances-list' in sys.argv or '-I' in sys.argv ))
    parser.add_argument('-I' ,'--instances-list', help='Path to the file with a list of the paths of the instances.', required=not ('-i' in sys.argv or '--instance' in sys.argv))
    parser.add_argument('-n' , default=1, type=int, help='Number of times each language will run opa opa opa')
    parser.add_argument('--lang' , nargs='+', required=True, help='Sources: python3, java, mcs, dotnet, julia, cpp, lua, javascript, matlab, golang')    
    parser.add_argument('--out' ,  default=data_dir,  help='Output dir')
    parser.add_argument('--timeout', default=DEFAULT_TIMEOUT, type=int,
            help='Timeout em segundos por execucao individual; ao estourar, pula as '
                 'repeticoes restantes e as instancias maiores dessa linguagem nesta sessao '
                 '(assume inst-list ordenada por tamanho crescente)')
    args = parser.parse_args()
    data_dir = args.out
    timeout = args.timeout

    sources = ["java", "dotnet", "mcs", "python3", "pypy", "julia", "cpp", "cppOOP",
            "fortran", "node", "lua", "luajit", "rust", "octave", "c", "matlab", "golang"]

    for i in args.lang:
        if i not in sources:
            print("{} is not suported".format(i))
            exit(0)

    sources = args.lang[:]
    instance = args.instance
    instances = []
    if args.instance != None:
        instances.append(args.instance)
    else:
        with open(args.instances_list, 'r') as f:
            instances = [line.strip() for line in f if line.strip()]

    n = args.n

    lang_dir = {
            "dotnet": "csharp",
            "mcs": "csharp",
            "java": "java",
            "python3": "python",
            "pypy": "python",
            "julia": "julia",
            "cpp" : "cplusplus",
            "cppOOP" : "cppOOP",
            "c" : "c",
            "fortran" : "fortran",
            "node" : "javascript",
            "lua" : "lua",
            "luajit" : "lua",
            "rust" : "rust",
            #"octave" : "octave",
            "matlab" : "octave",
            "golang" : "go"
            }

    """
    pnames = {
            "java" : "java",
            "csharp" : "cli",
            }
    """

    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    
    os.chdir("mlp-instances/loader")
    os.system("make")
    os.chdir("../../")

    timed_out_langs = set()

    for inst in instances:
        os.chdir("mlp-instances/")
        os.system("./load " + inst)
        os.chdir("../")

        for lang in sources:
            if lang in timed_out_langs:
                print("[skip-ahead] {} ja estourou timeout numa instancia menor, pulando {}".format(lang, inst))
                continue

            ds = ds_open(lang_dir[lang], data_dir)

            os.chdir(lang_dir[lang])
            print(os.getcwd())

            # IDEIA: criar arquivo cm lista de sources no dir de cada lang p gerar o nome do respectivo shell script "build_<source_name>.sh"
            if os.path.isfile("./build.sh"):
                # compila oq eh necessario
                os.system("./build.sh")

            for i in range(n):
                if lang == 'rust':
                    os.system("./build.sh")

                info = get_info(lang, timeout)
                info.update({"source" : lang, "instance" : inst, "branch" : get_branch()})

                ds = pd.concat([ds, pd.DataFrame(info)], ignore_index=True)

                print(ds)

                if info["timed_out"]:
                    timed_out_langs.add(lang)
                    print("[timeout] {} estourou {}s em {} - pulando repeticoes restantes e instancias maiores".format(lang, timeout, inst))
                    break

            ds.to_csv(os.path.join('..', data_dir,  lang_dir[lang] + '.csv'), index=False)

            os.chdir("..")

main()

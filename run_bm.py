import subprocess
import tempfile
import os
import pandas as pd
import concurrent.futures
import psutil as pU
import time

def ds_open(lang_name):
    ds_name = 'data/' + lang_name + ".csv"

    if os.path.isfile(ds_name):
        df = pd.read_csv(ds_name)
    else:
        os.system('cp data/template.csv data/' + lang_name + '.csv')
        df = pd.read_csv(ds_name)

    return df

def get_mem_avg(pid):

    p = pU.Process(pid=pid)
    interval = time.time()
    MB = 1024**2

    gap = 1e-10
    
    mem_max = -99999
    mem_avg = 0
    count = 0
    while True:
    #while count < 50:
        # time.sleep(0.1)
        #time.sleep(gap)

        # funciona durante 10s
        #print(time.time() - interval)
        #if time.time() - interval >= 10:
            #break

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

            """
            print(" Shared"  , float(mem.shared) / 1024**2)
            print(" Resident"  , float(mem.rss)/ 1024**2)
            print(" Virtual"  , float(mem.vms)/ 1024**2)
            print("pid", pid)
"""

            mem_total = (float(mem.rss)- float(mem.shared))/ MB
            if mem_total < gap:
                print("break by gap")
                break;
            elif mem_total > mem_max:
                mem_max = mem_total

            mem_avg += mem_total
        count += 1

    mem_avg /= count

    return {"mem_avg" : mem_avg, "mem_max" : mem_max, "mem_lookups" : count}

def get_COST(_str):
    lines = _str.split("\n")

    for l in lines :
        print(l)
        if l.find("COST") != -1:
            i = l.find(":")
            col = l[:i]
            val = float(l[i+1:].replace(',', '.'))
            return val

def get_TIME(_str):
    lines = _str.split("\n")

    for l in lines :
        print(l)
        if l.find("TIME") != -1:
            i = l.find(":")
            col = l[:i]
            val = float(l[i+1:].replace(',', '.'))
            return val
    return

def get_info(lang):
    f = open('run_' + lang + '.sh', 'r')
    print(f.readline())
    cmd_line = f.readline().replace('\n', '').split(' ')
    print(cmd_line)

    with tempfile.TemporaryFile() as tempf:
        proc = subprocess.Popen(cmd_line, stdout=tempf)
        pid = proc.pid

        info = get_mem_avg(pid)

        proc.wait()
        tempf.seek(0)
        output = str(tempf.read().decode("utf-8"))
        
        COST = get_COST(output)
        TIME = get_TIME(output)
        info.update({"COST" : [COST], "TIME" : [TIME]})
        print(info)

        return info

    f.close()

def main():

    # receber argumentos do usuario

    instance = "instances/att48.tsp"
    n = 2

    sources = ["julia"]#, "dotnet", "mcs", "python3", "pypy", "julia"]
    #sources = ["java", "dotnet", "mcs", "python3", "pypy", "julia"]
    lang_dir = {
            "dotnet": "csharp",
            "mcs": "csharp",
            "java": "java",
            "python3": "python",
            "pypy": "python",
            "julia": "julia"
            }

    """
    pnames = {
            "java" : "java",
            "csharp" : "cli",
            }
    """
    
    os.system("./load " + instance)

    for lang in sources:
        ds = ds_open(lang_dir[lang])
        
        os.chdir(lang_dir[lang])
        print(os.getcwd())

        # IDEIA: criar arquivo cm lista de sources no dir de cada lang p gerar o nome do respectivo shell script "build_<source_name>.sh"
        if os.path.isfile("./build.sh"):
            # compila oq eh necessario
            os.system("./build.sh")

        for i in range(n):

            info = get_info(lang)
            info.update({"source" : lang, "instance" : instance})

            ds = ds.append(pd.DataFrame(info), ignore_index=True)

            print(ds)

        ds.to_csv('../data/' +  lang_dir[lang] + '.csv', index=False)

        os.chdir("..")

main()

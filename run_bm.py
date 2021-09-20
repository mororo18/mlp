import subprocess
import tempfile
import os
import pandas as pd

def ds_open(lang_name):
    ds_name = 'data/' + lang_name + ".csv"

    try:
        df = pd.read_csv(ds_name)
    except FileNotFoundError:
        os.system('cp data/template.csv data/' + lang_name + '.csv')
        df = pd.read_csv(ds_name)

    return df

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

def get_info():
    with tempfile.TemporaryFile() as tempf:
        proc = subprocess.Popen(['./run.sh'], stdout=tempf)
        proc.wait()
        tempf.seek(0)
        output = str(tempf.read().decode("utf-8"))
        
        COST = get_COST(output)
        TIME = get_TIME(output)
        a = {"COST" : [COST], "TIME" : [TIME]}
        print(a)

        return a

def main():

    # receber argumentos do usuario

    n = 1

    langs_dirs = ["csharp"]
    #langs_dirs = ["java"]

    for lang in langs_dirs:
        ds = ds_open(lang)
        
        os.chdir(lang)
        print(os.getcwd())

        if os.path.isfile("./build.sh"):
            # compila oq eh necessario
            os.system("./build.sh")

        for i in range(n):
            """
            info = {
                "COST" :  ,
                "TIME" :  ,
                "mem"  :  ,
                "gpu"  :
            }
            """

            info = get_info()
            ds = ds.append(pd.DataFrame(info), ignore_index=True)
            print(ds)

        ds.to_csv('../data/' + lang + '.csv', index=False)

main()

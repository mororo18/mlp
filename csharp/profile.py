#! /usr/bin/python3

import subprocess
import signal
import os

#os.system('export DOTNET_PerfMapEnabled=1 ; export DOTNET_EnableEventLog=1')
#exit(0)
cs = subprocess.Popen(['dotnet','run'])

_id =  cs.pid

print(_id)

p = subprocess.Popen(['sudo', 'perf','record', '-p', str(_id), '-g', '-F' , '997'])

cs.wait()
p.send_signal(signal.SIGINT)


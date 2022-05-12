import os
import sys
import subprocess
import timeit
from pathlib import Path
import threading 
import queue
import itertools

cmds = [("Us","ptn.nwk","python3 /home/ewinship/project5.py "),("RapidNJ","rnj.nwk","/home/ewinship/rapidNJ/bin/rapidnj -i pd "),("Quicktree","qt.nwk","/home/ewinship/quicktree-2.0/quicktree -i m ")]

output = list()
my_q = queue.Queue()

class Worker(threading.Thread):
    def __init__(self,q):
        self.queue = q
        threading.Thread.__init__(self)
    
    def run(self):        
        while job := self.queue.get():            
            runinfo, file = job            
            name,sfx,cmd = runinfo
            if ".phy" not in file:
                self.queue.task_done()
                continue
            file_length = len(open(file).read().split("\n"))-1
            outfile = str(Path(file).with_suffix(''))
            outfile = outfile + "." + sfx
            name,sfx,cmd = runinfo 
            print(cmd+" "+file+" > "+outfile )
            timeit_str = "os.system('"+cmd+" "+file+" > "+outfile + "')" #os.system("python3 /home/ewinship/project5.py something.phy > something.ptn.nwk")
            ti = timeit.timeit(timeit_str,number=2,setup="import os")
            output.append((name, file_length,ti))
            self.queue.task_done()

jobs = list(itertools.product(cmds,os.listdir(sys.argv[1])))
#print(jobs)
for job in jobs:
    my_q.put(job)

threads = list()
for i in range(4):
    t = Worker(my_q)
    t.start()
    threads.append(t)

my_q.join()

for t in threads:
    threads.join()

with open("times2.txt","w") as f:
    for p,l,t in output:
        f.write("{}\t{:d}\t{:f}".format(p,l,t))

# for file in os.listdir(sys.argv[1]):
#     if ".phy" not in file:
#         continue
#     file_length = len(open(file).read().split("\n"))-1
#     for name,sfx,cmd in cmds:        
#         outfile = str(Path(file).with_suffix(''))
#         outfile = outfile + "." + sfx
#         print(cmd+" "+file+" > "+outfile )
#         timeit_str = "os.system('"+cmd+" "+file+" > "+outfile + "')"
#         ti = timeit.timeit(timeit_str,number=2,setup="import os")
#         programs.append(name)
#         lengths.append(file_length)
#         times.append(ti)
        


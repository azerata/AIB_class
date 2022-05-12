import timeit
import subprocess
import itertools
import os
import sys
import pandas as pd
import seaborn as sns
from pathlib import Path

###("DiscountNJ","dnj.nwk","python3 /home/hammer/project_05/discountNJ.py "),
cmds = [("RapidNJ","rnj.nwk","/home/hammer/project_05/rapidnj -i pd "),("Quicktree","qt.nwk","/home/hammer/project_05/quicktree -i m ")]

programs, lengths, times = [], [], []

#jobs = list(itertools.product(cmds,os.listdir(sys.argv[1])))
path = sys.argv[1]

for file in os.listdir(path):
    if ".phy" not in file:
        continue
    file_length = len(open(path+'/'+file).read().split("\n"))-1
    for name,sfx,cmd in cmds:        
        outfile = str(Path(file).with_suffix(''))
        outfile = outfile + "." + sfx
        print(cmd+" "+file+" > "+outfile )
        timeit_str = "os.system('"+cmd+" "+path+'/'+file+" > "+path+'/'+outfile + "')"
        ti = timeit.timeit(timeit_str,number=2,setup="import os")
        programs.append(name)
        lengths.append(file_length)
        times.append(ti)

time_measures = pd.DataFrame({'Program': programs, 'Length': lengths, 'time': times})
time_measures.to_csv("times.csv", sep = ",")

g = sns.lmplot(x = 'Length', y = 'time', hue = 'Program', lowess = True,
           x_jitter = 0.1, markers='.',
           data = time_measures)
g.set(ylim = (0, max(time_measures['time'])))

g.savefig(path+"/"+"comparison.png")

dist_list = []
for file in os.listdir(path):
    if ".nwk" not in file:
        continue
    print(file)
    dist_list.append(file)
print(dist_list)

file_1, file_2, rfdist = [], [], []

leng = len(dist_list)-1
i = leng
while i >= 0:
    j = leng
    while j >= 0:
        if i == j or i < j:
            pass
        else:
            f_1 = dist_list[i]
            f_2 = dist_list[j]
            if f_1[0:2] == f_2[0:2]:
                file_1.append(dist_list[i])
                file_2.append(dist_list[j])
                dist = subprocess.run(["python3", "rfdist_1.py", f"{path}/{f_1}", f"{path}/{f_2}"], capture_output=True,text=True)
                print(dist.stderr)
                rfdist.append(dist.stdout)
        j-=1
    i-=1


with open(path+"/"+"times2.txt","w") as f:
    for p in range(leng):
        f.write("{}\t{}\t{}\n".format(file_1[p],file_2[p],rfdist[p]))
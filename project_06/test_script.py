from re import sub
import subprocess
import sys
import timeit
import pandas as pd
import seaborn as sns




num, n, s, optimal = [],[],[],[] 

with open(sys.argv[1]) as file:
    for line in file.readlines():
        tmp = line.split()
        num.append(tmp[0])
        n.append(len(tmp[1]))
        s.append(tmp[1])
        optimal.append(tmp[2])

df = pd.DataFrame({"seq #": num, "length": n, "sequence": s, "Optimal score": optimal})


time, fold = [], []

for i in enumerate(df['sequence']):
    cmd = f"subprocess.run(['python3', 'hp1_4.py', '{i[1]}'], capture_output=True,text=True)"
    ti = timeit.timeit(cmd,number=2,setup="import subprocess")
    outp = subprocess.run(['python3', 'hp1_4.py', f'{i[1]}'], capture_output=True,text=True)
    fold.append(outp.stdout.strip('\n'))
    time.append(ti)

tmp = pd.DataFrame({"Times":time, "Fold":fold})

df2 = pd.concat([df,tmp], axis=1)


#print(df2)

ds = []

for i in enumerate(df2["sequence"]):
    tmp = subprocess.run(["python3","hpview3k.py", f'{df2["sequence"][i[0]]}', f'{df2["Fold"][i[0]]}'],capture_output=True, text=True)
    score = tmp.stdout[-4::].strip(':')
    ds.append(score.strip('\n'))

df3 = pd.DataFrame({"Discovered score":ds})
final = pd.concat([df2, df3],axis=1)

final.to_csv("output.csv", sep = ",")
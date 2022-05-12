import typing

def read_string(path:str)->str:
    with open(path) as file
    out = file.readline()
    return out

def match_oe(x:str)-> list:
    match_1, match_2 = [], []
    score_1, score_2 = 0, 0

    i = 0
    j = len(x)
    while i < j:
        if i%2==j%2:
            pass
        elif i % 2 == 0 and j % 2 != 0:
            if x[i] in ["h", "H"] and x[j] in ["h", "H"]:
                match_1.append((i,j))
                score_1+=1
        else: 
            if x[i] in ["h", "H"] and x[j] in ["h", "H"]:
                match_2.append((i,j))
                score_2+=1
        i+=1
        j-=1
    


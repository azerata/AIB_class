import typing
import sys

def read_string(path:str)->str:
    with open(path) as file:
        out = [line for line in file.readlines()]
        out = ''.join(out)
        return out

def match_oe(x:str)-> list:
    match_1, match_2 = [], []
    score_1, score_2 = 0, 0
    i = 0
    j = len(x)-1
    while i < j:
        if i % 2 == 0:
            if j % 2 != 0:
                if x[i] in ['h', 'H']:
                    if x[j] in ['h', 'H']:
                        match_1.append((i,j))
                        score_1-=1
                        i+=2
                        j-=2
                    else:
                        j-=2
                else:
                    i+=2
            else:
                j-=1
        else:
            i+=1
        
    i = 1
    j = len(x)-1
    while i < j:
        if i % 2 != 0:
            if j % 2 == 0:
                if x[i] in ['h', 'H']:
                    if x[j] in ['h', 'H']:
                        match_2.append((i,j))
                        score_2-=1
                        i+=2
                        j-=2
                    else:
                        j-=2
                else:
                    i+=2
            else:
                j-=1
        else:
            i+=1


    return [(score_1, match_1), (score_2, match_2)]

def merger(l1:list, l2:list):
    out = []
    
    i=j=0
    while i < len(l1) and j < len(l2):
        if l1[i][0] < l2[j][0]:
            out.append(l1[i])
            i+=1
        else:
            out.append(l2[j])
            j+=1
    if i < len(l1):
        out.extend(l1[i:])
    else:
        out.extend(l2[j:])
    return out

def calculate_u(direction, dist):
    if direction == 'e':
        if dist % 2 == 0:
            height = (dist // 2)-1
            out:list = height * ['n'] + ['e'] + height * ['s'] +['e']
        if dist % 2 != 0:
            height = (dist // 2)
            out:list = height * ['n'] + ['e'] + height * ['s'] 
    if direction == 'w':
        if dist % 2 == 0:
            height = (dist // 2)-1
            out:list = height * ['s'] + ['w'] + height * ['n'] +['w']
        if dist % 2 != 0:
            height = (dist // 2)
            out:list = height * ['s'] + ['w'] + height * ['n'] 
    if direction == 't':
        if dist % 2 == 0:
            height = (dist // 2)-1
            out:list = height * ['e'] + ['s'] + height * ['w'] +['w']
        if dist % 2 != 0:
            height = (dist // 2)
            out:list = height * ['e'] + ['s'] + height * ['w']

    return out


def calculate_fold(x, seq):
    out = []
    matches = len(x)
    for i in range(matches):
        if i < matches -1:
            distance = x[i+1][0] - x[i][0]
            if distance < 3:
                tmp = distance * ["e"]
            else:
                tmp = calculate_u('e', distance)
        else:
            distance = x[i][1] - x[i][0]
            tmp = calculate_u('t', distance)
        out.extend(tmp)
    for i in range(matches-1, 0,-1):
        if i > 0:
            distance = x[i-1][1] - x[i][1]
            if distance < 3:
                tmp = distance * ['w']
            else:
                tmp = calculate_u('w', distance)
        out.extend(tmp)
    if len(seq) > len(out):
        tmp = (len(seq)-x[0][1]-1) * ['w']
        out.extend(tmp)
        tmp = x[0][0] * ['e']
        tmp.extend(out)
        out = tmp
    return out

def fold_mark_II(x, seq):
    out = []
    matches = len(x)
    start = x[-1]
    out.extend(calculate_u('t',start[1] - start[0] ))
    n = []
    s = []
    for i in enumerate(x):
        if i[0] < matches-1:
            n.extend(calculate_u('e', ))
    return out
    

    
    

if __name__=='__main__':
    seq = read_string(sys.argv[1])
    scores = match_oe(seq)
    

    fold = calculate_fold(scores[0][1], seq)
    print(''.join(fold))
    print(len(seq), len(fold))
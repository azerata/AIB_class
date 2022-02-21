import numpy as np
import sys
from cmath import inf


tmp = sys.stdin.readline()
score = np.zeros(tmp,tmp)
scores = [sys.stdin.readline() for i in range(tmp)]

print(scores)

def cost(i:str, j:str, score = None )->int:
    ''' 
    Score matrix for match/mismatch
    With translation from Nucleotide to idx
    Takes two nucleotides as input (strings of length 1)
    And score matrix, (numpy array)
    outputs the score for the given match/mismatch as int
    '''


    t = {"A":0,
         "C":1,
         "G":2,
         "T":3}
    return score[t[i],t[j]]
#%%
from cmath import inf
import numpy
import numpy as np
from fasta import fasta

'''
Dependant on packages: numpy & math
Also needs CHR's fasta parser in same dir as this is run from 
That and the input files ofcourse... 
'''


def cost(i:str, j:str)->int:
    ''' 
    Score matrix for match/mismatch
    With translation from Nucleotide to idx
    Takes to nucleotides as input (strings of length 1)
    outputs the score for the given match/mismatch as int
    '''

    score = numpy.array([[0, 5, 2, 5],
                         [5, 0, 5, 2],
                         [2, 5, 0, 5],
                         [5, 2, 5, 0]])

    t = {"a":0,
         "c":1,
         "g":2,
         "t":3}
    return score[t[i],t[j]]

#%%

#%% build matrix

def AffineScore(seq1,seq2,a,b):
    S=np.zeros((len(seq1)+1,len(seq2)+1))
    D=np.zeros((len(seq1)+1,len(seq2)+1))
    I=np.zeros((len(seq1)+1,len(seq2)+1))

    for i in range(0,len(seq1)+1):
        for j in range (0,len(seq2)+1):

            #Compute D[i,j]
            v1,v2 = inf, inf
            if i > 0 and j>=0:
                v1=S[i-1,j]+(a+b)
            if i > 1 and j >=0:
                v2=D[i-1,j]+a
            D[i,j]=min(v1,v2)

            #Compute I[i,j]
            v1 = v2 = inf
            if i >= 0 and j > 0:
                v1=S[i,j-1]+(a+b)
            if i >= 0 and j > 1:
                v2=I[i,j-1]+a
            I[i,j]=min(v1,v2)

            #Compute  S[i,j]
            v1 = v2 =v3 = v4 = inf
            if i == 0 and j == 0:
                v1 = 0
            if i > 0 and j > 0:
                v2 = S[i-1,j-1]+cost(seq1[i-1],seq2[j-1])
            if i > 0 and j >= 0:
                v3 = D[i,j]
            if i >= 0 and j > 0:
                v4 = I[i,j]
            S[i,j] = min(v1,v2,v3,v4)

    return [S,D,I]

def global_affine (seq1,seq2,a,b):
    S = AffineScore(seq1,seq2,a,b)[0]
    str1 = []
    str2 = []
    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i,j] == S[i-1,j-1]+cost(seq1[i-1],seq2[j-1]):
            str1.append(seq1[i-1])
            str2.append(seq2[j-1])
            i=i-1
            j=j-1
        else:
            k=1
            while True:
                if i >= k and S[i,j]==S[i-k,j]+(a*k+b):
                    str1.append(seq1[i-1-k:i-1][::-1])
                    str2.append("-"*k)
                    i=i-k
                    break

                elif j >k and S[i,j] == S[i,j-k]+(a*k+b):
                    str1.append("-"*k)
                    str2.append(seq2[j-1-k:j-1][::-1])
                    j=j-k
                    break

                else:
                    k=k+1
    return "".join(str1[::-1]),"".join(str2[::-1])




'''
Area below for testing:
'''

#%% Test block 1 (Task 1)

#Case1
#seq2 = "ACGTCGTAGCTA"
#seq1 = "ACGTGTCAACGT"

#Case2
#seq2 = "AAGG"
#seq1 = "AATAAT"

#Case3
#seq1 = "TCCAGAGA"
#seq2 = "TCGAT"

#Case4
seq1 = "ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtctgtttacgtataaacagaatcgcctgggttcgc"
seq2 = "gggctaaaggttagggtctttcacactaaagagtggtgcgtatcgtggctaatgtaccgcttctggtatcgtggcttacggccagacctacaagtactagacctgagaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc"


print(AffineScore(seq1,seq2,5,5))
print(global_affine (seq1,seq2,5,5))




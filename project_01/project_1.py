#%%
from cmath import inf
import numpy
from fasta import fasta
gap_cost = -5

def cost(i:str, j:str)->int:
    score = numpy.array([[10, 2, 5, 2],
                         [2, 10, 2, 5],
                         [5, 2, 10, 2],
                         [2, 5, 2, 10]])

    t = {"A":0,
         "C":1,
         "G":2,
         "T":3}
    return score[t[i],t[j]]

    
#%%

#%% build matrix
seq2 =  "AATAAT"
seq1 = "AAGG"


i,j = len(seq1),len(seq2)

def setup_matrix(seq1, seq2):
    matrix = numpy.full([len(seq1)+1,len(seq2)+1], None)
    return matrix


# %%
def c(i, j, matrix, gap_cost):
    if matrix[i,j] is None:
        v1, v2, v3, v4 = -inf, -inf, -inf, -inf
        if i > 0 and j > 0:
            v1 = c(i-1,j-1, matrix, gap_cost) + cost(seq1[i-1],seq2[j-1])
        if i > 0 and j >= 0:
            v2 = c(i-1,j, matrix, gap_cost) + gap_cost
        if i >= 0 and j > 0:
            v3 = c(i, j-1, matrix, gap_cost) + gap_cost
        if i == 0 and j == 0:
            v4 = 0
        matrix[i,j] = max(v1, v2, v3, v4)
    return matrix[i,j]

#%%
def build_alignment(seq1, seq2, gap_score):
    matrix = setup_matrix(seq1, seq2)
    i, j = len(seq1), len(seq2)
    out = c(i, j, matrix, gap_score)
    return out, matrix


#%%
print(build_alignment(seq1, seq2, -5))

# %%

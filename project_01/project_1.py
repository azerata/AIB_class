#%%
from cmath import inf
import numpy

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

matrix = numpy.full([len(seq1)+1,len(seq2)+1], None)

i,j = len(seq1),len(seq2)



# %%
def c(i, j):
    if matrix[i,j] is None:
        v1, v2, v3, v4 = -inf, -inf, -inf, -inf
        if i > 0 and j > 0:
            v1 = c(i-1,j-1) + cost(seq1[i-1],seq2[j-1])
        if i > 0 and j >= 0:
            v2 = c(i-1,j) + gap_cost
        if i >= 0 and j > 0:
            v3 = c(i, j-1) + gap_cost
        if i == 0 and j == 0:
            v4 = 0
        matrix[i,j] = max(v1, v2, v3, v4)
        print(matrix)
    return matrix[i,j]

#%%


#%%
mat = c(i, j)
print(matrix)

# %%

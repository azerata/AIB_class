#%%
from cmath import inf
import numpy
import sys
from fasta import fasta
'''
Dependant on packages: numpy & math
Also needs CHR's fasta parser in same dir as this is run from 
That and the input files ofcourse... 
'''


def cost(i:str, j:str, score = None )->int:
    ''' 
    Score matrix for match/mismatch
    With translation from Nucleotide to idx
    Takes two nucleotides as input (strings of length 1)
    And score matrix, (numpy array)
    outputs the score for the given match/mismatch as int
    '''

    if score is None:
        score = numpy.array([[10, 2, 5, 2],
                             [2, 10, 2, 5],
                             [5, 2, 10, 2],
                             [2, 5, 2, 10]])

    t = {"A":0,
         "C":1,
         "G":2,
         "T":3}
    return score[t[i],t[j]]


#%% build matrix


def setup_matrix(seq1:str, seq2:str):
    '''
    Takes two strings as input
    Strings should be two nucleotide sequences we wish to align
    Returns a matrix of size seq1 X seq2

    '''

    matrix = numpy.full([len(seq1)+1,len(seq2)+1], None)
    return matrix


# %%
def c(i:int, j:int, matrix, gap_cost:int, seq1:str, seq2:str, score)->int:
    ''' 
    Input: 
    Two integers i and j, each equivalent to the length of the corresponding nucleotide sequence
    A matrix to work on (numpy array)  as output by setup_matrix function
    A linear gap modifier (int) and 2 sequences

    Returns the optimal alignment cost of the two sequenses
    '''

    if matrix[i,j] is None:
        v1, v2, v3, v4 = -inf, -inf, -inf, -inf
        if i > 0 and j > 0:
            v1 = c(i-1,j-1, matrix, gap_cost, seq1, seq2, score) + cost(seq1[i-1],seq2[j-1], score)
        if i > 0 and j >= 0:
            v2 = c(i-1,j, matrix, gap_cost, seq1, seq2, score) + gap_cost
        if i >= 0 and j > 0:
            v3 = c(i, j-1, matrix, gap_cost, seq1, seq2, score) + gap_cost
        if i == 0 and j == 0:
            v4 = 0
        matrix[i,j] = max(v1, v2, v3, v4)
    return matrix[i,j]

#%%
def build_alignment(seq1:str, seq2:str, gap_score:int, score = None):
    '''
    Wrapper for the previous functions, call them in order to pipe inputs/outputs
    Takes 3 inputs, 2 sequences, and a gap score
    returns a tuple with the filled out matrix, and the an int the optimal alignment cost
    '''

    matrix = setup_matrix(seq1, seq2)
    i, j = len(seq1), len(seq2)
    out = c(i, j, matrix, gap_score, seq1, seq2, score)
    return  (matrix, out)




if __name__ == '__main__':
    
    seqs = sys.stdin.readlines()
    build_alignment(seqs[0],seqs[1], -5)



'''
Area below for testing:
'''

#%% Test block 1 (Task 1)

seq2 =  "AATAAT"
seq1 = "AAGG"

build_alignment(seq1, seq2, -5)

# %% Test block 2 (fasta seqs, task 2)
o = fasta("seq1.fasta")
o2 = fasta("seq2.fasta")

seq1 = o["Seq1"]
seq2 = o2["Seq2"]
i = build_alignment(seq1, seq2, -5)
print(i[0])
print(i[1])

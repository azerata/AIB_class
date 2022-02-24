import numpy as np
import sys
import pandas as pd
from timeit import default_timer as timer
import seaborn as sns

from fasta import fasta
from cmath import inf



def outer_cost():

    file = open(sys.argv[1])
    tmp = int(file.readline())

    t = {}
    score = np.zeros((tmp,tmp))
    scores = [file.readline() for i in range(tmp)]

    i = 0
    for val in scores:
        tmp = [_ for _ in val.split()]
        t[tmp[0]]=i
        for j in range(len(tmp[1:])):
            score[i,j] = int(tmp[j +1])
        i+=1
    file.close()
    def cost(i,j): return score[t[i],t[j]]

    return cost

def setup_matrix(seq1:str, seq2:str):
    '''
    Takes two strings as input
    Strings should be two nucleotide sequences we wish to align
    Returns a matrix of size seq1 X seq2

    '''

    matrix = np.full([len(seq1)+1,len(seq2)+1], None)
    return matrix

def c(i:int, j:int, matrix, gap_cost:int, seq1:str, seq2:str)->int:
    ''' 
    Input: 
    Two integers i and j, each equivalent to the length of the corresponding nucleotide sequence
    A matrix to work on (numpy array)  as output by setup_matrix function
    A linear gap modifier (int) and 2 sequences

    Returns the optimal alignment cost of the two sequenses
    '''

    if matrix[i,j] is None:
        v1, v2, v3, v4 = inf, inf, inf, inf
        if i > 0 and j > 0:
            v1 = c(i-1,j-1, matrix, gap_cost, seq1, seq2) + cost(seq1[i-1],seq2[j-1])
        if i > 0 and j >= 0:
            v2 = c(i-1,j, matrix, gap_cost, seq1, seq2) + gap_cost
        if i >= 0 and j > 0:
            v3 = c(i, j-1, matrix, gap_cost, seq1, seq2) + gap_cost
        if i == 0 and j == 0:
            v4 = 0
        matrix[i,j] = min(v1, v2, v3, v4)
    return matrix[i,j]

def build_alignment(seq1:str, seq2:str, gap_score:int):
    '''
    Wrapper for the previous functions, call them in order to pipe inputs/outputs
    Takes 3 inputs, 2 sequences, and a gap score
    returns a tuple with the filled out matrix, and the an int the optimal alignment cost
    '''

    matrix = setup_matrix(seq1, seq2)
    i, j = len(seq1), len(seq2)
    out = c(i, j, matrix, gap_score, seq1, seq2)
    return  (matrix, out)

def backtrack(matrix, i,j, seq1, seq2, gap_cost, out1 = [], out2 = []):
    T = matrix
    if (i > 0) and (j > 0) and (T[i,j] == (T[i-1,j-1] + cost(seq1[i-1],seq2[j-1]))):
        out1.append(seq1[i-1])
        out2.append(seq2[j-1])
        backtrack(T, i-1, j-1, seq1, seq2, gap_cost, out1, out2)
    elif (i > 0) and (j >= 0) and (T[i,j] == (T[i-1,j] + gap_cost)):
        out1.append(seq1[i-1])
        out2.append("-")
        backtrack(T, i-1, j, seq1, seq2, gap_cost, out1, out2)
    elif (i >= 0) and (j > 0) and (T[i,j] == (T[i,j-1] + gap_cost)):
        out1.append("-")
        out2.append(seq2[j-1])
        backtrack(T, i, j-1, seq1, seq2, gap_cost, out1, out2)
    return("".join(out1)[::-1], "".join(out2)[::-1] )

def optimal_alignment(seq1, seq2, gap_cost):
    pip = build_alignment(seq1, seq2, gap_cost)
    out = backtrack(pip[0], len(seq1), len(seq2), seq1, seq2, gap_cost)
    return out



def get_sequences(file):
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:
        out.append(seq_dict[key].upper())
    return out

def optimal_score_matrix(sequences:list, gap_cost):
    k = len(sequences)
    M = np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            S = build_alignment(sequences[i],sequences[j],gap_cost)[1]
            M[i,j] = S
    return M


def construct_test_data(sequences:list, n):

    out = [ [] for _ in range(n)]

    for i in range(n):
        for j in range(len(sequences)):
            tmp = len(sequences[j]) // n
            out[i].append(sequences[j][:(i+1)*tmp])

    return out



def time_alg(test_data):

    ns, times = [], []
    

    for n in range( len(test_data)):
        start = timer()
        optimal_alignment(test_data[n][0], test_data[n][1],gap_cost )
        end = timer()

        ns.append(len(test_data[n][0]))
        times.append(end - start)

    return pd.DataFrame({'n':ns, 'time':times})

if __name__ == "__main__":

    sequences = get_sequences(sys.argv[2])
    k = len(sequences)

    cost = outer_cost()

    gap_cost = int(sys.argv[3])

    test_data = construct_test_data(sequences,20)

    out = time_alg(test_data)
    out.to_excel("test_1.xlsx")

#    if k == 2:
#        out = optimal_alignment(sequences[0], sequences[1], gap_cost)
#        print(out[0], "\n", out[1])
#    elif k > 2:
#        out = optimal_score_matrix(sequences, gap_cost)
#        print(out)


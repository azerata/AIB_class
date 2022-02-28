from cgi import test
import numpy as np
import pandas as pd
from timeit import default_timer as timer
import seaborn as sns
import sys

from cmath import inf
from fasta import fasta

def construct_test_data(sequences:list, n):

    out = [ [] for _ in range(n)]

    for i in range(n):
        for j in range(len(sequences)):
            tmp = len(sequences[j]) // n
            out[i].append(sequences[j][:(i+1)*tmp])
    f = [len(seq[0]) for seq in out]
    print(f"test data constructed, lengths: {f}")
    return out


def time_affine(test_data):

    ns, times, alg = [], [], []
    

    for n in range(len(test_data)):
        start = timer()
        global_affine(test_data[n][0], test_data[n][1],gap_cost["alpha"], gap_cost["beta"] )
        end = timer()
        print(f"Completed run: {n+1}/{len(test_data)}, at seq lengths {len(test_data[n][0])} & {len(test_data[n][1])}.")

        ns.append(len(test_data[n][0]))
        times.append(end - start)
        alg.append("affine")

    return pd.DataFrame({'n':ns, 'time':times,'algorithm':alg})


def time_linear(test_data):

    ns, times, alg = [], [], []
    

    for n in range( len(test_data)):
        start = timer()
        optimal_alignment(test_data[n][0], test_data[n][1],gap_cost["alpha"] )
        end = timer()

        print(f"Completed run: {n+1}/{len(test_data)}, at seq lengths {len(test_data[n][0])} & {len(test_data[n][1])}.")
        ns.append(len(test_data[n][0]))
        times.append(end - start)
        alg.append("linear")

    return pd.DataFrame({'n':ns, 'time':times, 'algorithm':alg})


def outer_cost():
    '''
    Outer function to create the cost function from phylib format score matrix
    Builds a score matrics using NP arrays
    And a translation matrix
    Uses sys arg 1 as input file
    '''
    file = open(sys.argv[1])
    tmp = int(file.readline())

    t = {}
    score = np.zeros((tmp,tmp))
    scores = [file.readline() for i in range(tmp)]

    i = 0
    for val in scores:
        tmp = [_ for _ in val.split()]
        t[tmp[0].upper()]=i
        t[tmp[0].lower()]=i
        for j in range(len(tmp[1:])):
            score[i,j] = int(tmp[j +1])
        i+=1
    file.close()
    def cost(i,j): return score[t[i],t[j]]

    return cost


def get_sequences(file):
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:
        out.append(seq_dict[key].upper())
    return out

#### Linear ####

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


def linear_optimal_score_matrix(sequences:list, gap_cost):
    k = len(sequences)
    M = np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            S = build_alignment(sequences[i],sequences[j],gap_cost)[1]
            M[i,j] = S
    return M



### Affine ####

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


def affine_optimal_score_matrix(sequences:list, alpha:int, beta:int):
    k = len(sequences)
    M = np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            S = AffineScore(sequences[i], sequences[j], alpha, beta)[0]
            M[i,j] = S[len(sequences[i]), len(sequences[j])]
    return M







cost = outer_cost()

sequences = get_sequences(sys.argv[2])
sequences.sort(key=len,reverse=True)

gap_cost = {"alpha":int(sys.argv[3]), "beta":int(sys.argv[4])}


test_data = construct_test_data(sequences, 20)

t1 = time_affine(test_data)
t2 = time_linear(test_data)


time_measures = pd.concat([t1, t2])

g = sns.lmplot(x = 'n', y = 'time', hue='algorithm', lowess = True,
           x_jitter = 0.1, markers='.',
           data = time_measures)
g.set(ylim = (0, max(time_measures['time'])))

g.savefig("comparison.png")
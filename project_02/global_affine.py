from cmath import inf


import numpy as np
import sys

from fasta import fasta


def outer_cost():
    '''
    Outer function to create the cost function from phylib format score matrix
    Builds a score matrics using NP arrays
    And a translation matrix
    Uses sys arg 2 as input file
    '''
    file = open(sys.argv[2])
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


def optimal_score_matrix(sequences:list, alpha:int, beta:int):
    k = len(sequences)
    M = np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            S = AffineScore(sequences[i], sequences[j], alpha, beta)[0]
            M[i,j] = S[len(sequences[i]), len(sequences[j])]
    return M



def get_sequences(file):
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:
        out.append(seq_dict[key].upper())
    return out


if __name__ == "__main__":

    run_opt = sys.argv[1]

    cost = outer_cost() #uses sys.argv[2]

    sequences = get_sequences(sys.argv[3])

    gap = {"alpha":int(sys.argv[4]), "beta":int(sys.argv[5])}

    if run_opt in ['a', 'A']:
        print(global_affine(sequences[0], sequences[1], gap["alpha"], gap["beta"]))

    if run_opt in ['m', 'M']:
        print(optimal_score_matrix(sequences, gap["alpha"], gap["beta"]))
    
    if run_opt in ['o', 'O']:
        out = AffineScore(sequences[0], sequences[1], gap["alpha"], gap["beta"])[0][-1,-1]
        print(out)
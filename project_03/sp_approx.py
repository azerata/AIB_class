import numpy as np
import sys

from global_linear import optimal_alignment, optimal_score_matrix
from fasta import fasta


def outer_cost(path):
    '''
    Outer function to create the cost function from phylib format score matrix
    Builds a score matrics using NP arrays
    And a translation matrix
    
    Takes 1 input filepath, outputs a function to calculate cost. 
    '''
    file = open(path)
    tmp = int(file.readline())

    t = {"n":0, "N":0, "r":0, "R":0, "s":0, "S":0}
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

def get_sequences(file, lim = 'a'):
    ''' 
    Takes variable "file" that is the path to the wanted fasta file
    uses fasta function by storm, to parse, and outputs dict, in which all sequences can be found.
    '''
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:    
        out.append(seq_dict[key].upper())
    if lim not in ['a', 'A']:
        out = out[0:int(lim)]
    return out, keys


if __name__ == "__main__":

    cost = outer_cost(sys.argv[1])

    sequences, seq_names = get_sequences(sys.argv[2], sys.argv[4])

    gap = int(sys.argv[3])

    sys.setrecursionlimit(10**5)

    T = optimal_score_matrix(sequences, gap, cost)
    

    tmp = [ sum(row) for row in T ]
    
    s = np.argmin(tmp)

    k = len(sequences)

    i = 0

    MA = []
    while i < k:
        M = MA
        if i != s:
            seqs = optimal_alignment(sequences[s], sequences[i], gap, cost)
            A = []
            for n in zip(seqs[0], seqs[1]):
                A.append([n[0], n[1]])

            if  not MA:
                MA = A
            else:
                MA = []
                p = q = 0
                while p < len(M) and q < len(A):
                    if M[p][0] == '-' and A[q][0] == '-':
                        M[p].append(A[q][1])
                        MA.append(M[p])
                        p += 1
                        q += 1
                    elif M[p][0] == '-' and A[q][0] != '-':
                        M[p].append('-')
                        MA.append(M[p])
                        p += 1
                    elif M[p][0] != '-' and A[q][0] == '-':
                        tmp = ['-' for _ in range(len(M[p]))]
                        tmp.append(A[q][1])
                        MA.append(tmp)
                        q += 1
                    elif M[p][0] != '-' and A[q][0] != '-':
                        M[p].append(A[q][1])
                        MA.append(M[p])
                        p += 1
                        q += 1
                
                if p < len(M):
                    while p < len(M):
                        M[p].append('-') 
                        MA.append(M[p])
                        p += 1
                #elif q < len(A):
                #    while q < len(A):
                #        tmp = ['-' for _ in range(len(M[0]))]
                #        tmp[0] = A[q][0]
                #        tmp[-1] = A[q][1]
                #        MA.append(tmp)
                #        q += 1
        i +=1
        

    out =[ [] for _ in range(len(MA[0]))]

    for i in range(len(out)):
        for j in range(len(MA)):
            out[i].append(MA[j][i])


    printer = [ ''.join(seq) for seq in out]

    seq_list = []
    seq_list.append(seq_names[s])

    for i in range(0, k):
        if i != s:
            seq_list.append(seq_names[i])

    for n in zip(seq_list, printer):
        print(f'>{n[0]}',f'{n[1].strip()}',sep='\n')










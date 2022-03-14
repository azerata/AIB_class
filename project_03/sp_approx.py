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

def get_sequences(file):
    ''' 
    Takes variable "file" that is the path to the wanted fasta file
    uses fasta function by storm, to parse, and outputs dict, in which all sequences can be found.
    '''
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:
        out.append(seq_dict[key].upper())
    return out


if __name__ == "__main__":

    cost = outer_cost(sys.argv[1])

    sequences = get_sequences(sys.argv[2])

    gap = int(sys.argv[3])

    T = optimal_score_matrix(sequences, gap, cost)
    

    tmp = [sum(row) for row in T ]
    
    s = np.argmin(tmp)

    k = len(sequences)

    M = []

    i = 0

    seqs = optimal_alignment(sequences[s], sequences[i], gap, cost)
    A = []
    for n in zip(seqs[0], seqs[1]):
        A.append([n[0], n[1]])

    print(A)





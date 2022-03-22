from cmath import inf
import numpy
import numpy as np
import sys
from fasta import fasta


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
         "t":3,
         "A":0,
         "C":1,
         "G":2,
         "T":3}
    return score[t[i],t[j]]


def MSA_3(seq1,seq2,seq3,gap): 
	''' 
    Computing an optimal MSA using a column based score
    Iterative formulation for 3 seqs
    Takes 3 seqs as inputs and the gapcost
    Output the 3D score table
    '''
	T=np.zeros((len(seq1)+1,len(seq2)+1,len(seq3)+1))
	for i in range(len(seq1)+1):
		for j in range(len(seq2)+1):
			for k in range(len(seq3)+1):
				v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = inf
				if i==0 and j==0 and k==0:
					v0=0
				if i>0 and j>0 and k>0:
					v1=T[i-1][j-1][k-1]+cost(seq1[i-1],seq2[j-1])+cost(seq3[k-1],seq2[j-1])+cost(seq1[i-1],seq3[k-1])
				if i>0 and j>0 and k>=0: 
					v2=T[i-1][j-1][k] + cost(seq1[i-1],seq2[j-1])+gap+gap
				if i>0 and j>=0 and k>0:
					v3=T[i-1][j][k-1] + gap + cost(seq1[i-1],seq3[k-1])+gap
				if i>=0 and j>0 and k>0:
					v4 = T[i][j-1][k-1] + gap + gap + cost(seq2[j-1],seq3[k-1])
				if i>0 and j>=0 and k>=0:
					v5 = T[i-1][j][k] + gap + gap
				if i>=0 and j>0 and k>=0:
					v6 = T[i][j-1][k] + gap + gap
				if i>=0 and j>=0 and k>0:
					v7 = T[i][j][k-1] + gap + gap
				T[i][j][k]=min(v0,v1,v2,v3,v4,v5,v6,v7)

	return T




def backtrack(matrix, i,j,k, seq1, seq2, seq3, gap_cost, out1 = [], out2 = [], out3 = []):
	T = matrix

	if i>0 and j>0 and k>0 and (T[i,j,k]==T[i-1][j-1][k-1]+cost(seq1[i-1],seq2[j-1])+cost(seq3[k-1],seq2[j-1])+cost(seq1[i-1],seq3[k-1])):
		out1.append(seq1[i-1])
		out2.append(seq2[j-1])
		out3.append(seq3[k-1])
		backtrack(T, i-1, j-1,k-1, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	elif i>0 and j>0 and k>=0 and (T[i,j,k]==T[i-1][j-1][k] + cost(seq1[i-1],seq2[j-1])+gap+gap):
		out1.append(seq1[i-1])
		out2.append(seq2[j-1])
		out3.append("-")
		backtrack(T, i-1, j-1,k, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	elif i>0 and j>=0 and k>0 and (T[i,j,k]==T[i-1][j][k-1] + gap + cost(seq1[i-1],seq3[k-1])+gap):
		out1.append(seq1[i-1])
		out2.append("-")
		out3.append(seq3[k-1])
		backtrack(T, i-1, j,k-1, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	elif i>=0 and j>0 and k>0 and (T[i,j,k]==T[i][j-1][k-1] + gap + gap + cost(seq2[j-1],seq3[k-1])):
		out1.append("-")
		out2.append(seq2[j-1])
		out3.append(seq3[k-1])
		backtrack(T, i, j-1,k-1, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	elif i>0 and j>=0 and k>=0 and (T[i,j,k]== T[i-1][j][k] + gap + gap):
		out1.append(seq1[i-1])
		out2.append("-")
		out3.append("-")
		backtrack(T, i-1, j,k, seq1, seq2, seq3, gap_cost, out1, out2, out3)   
	elif i>=0 and j>0 and k>=0 and (T[i,j,k]==T[i][j-1][k] + gap + gap):
		out1.append("-")
		out2.append(seq2[j-1])
		out3.append("-")
		backtrack(T, i, j-1,k, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	elif i>=0 and j>=0 and k>0 and (T[i,j,k]==T[i][j][k-1] + gap + gap):
		out1.append("-")
		out2.append("-")
		out3.append(seq3[k-1])
		backtrack(T, i, j,k-1, seq1, seq2, seq3, gap_cost, out1, out2, out3)
	return("".join(out1)[::-1], "".join(out2)[::-1] ,"".join(out3)[::-1] )

def get_sequences(file):
    seq_dict = fasta(file)
    keys = [key for key in seq_dict.keys()]
    out = []
    for key in keys:
        out.append(seq_dict[key].upper())
    return out

if __name__ == "__main__":

	gap=5

	x=(get_sequences(sys.argv[1]))

	matrix=MSA_3(x[0],x[1],x[2],gap)

	print(matrix[-1][-1][-1])

	print(backtrack(matrix, len(x[0]),len(x[1]),len(x[2]), x[0], x[1], x[2], gap))


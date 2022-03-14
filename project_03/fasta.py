#
# fasta.py
#
# Christian Storm Pedersen

import sys
import string

def fasta(f):
    """
    Reads the fasta file f and returns a dictionary with the sequence names as keys and the
    sequences as the corresponding values. Lines starting with ';' in the fasta file are
    considered comments and ignored.
    """
    d = {}
    curr_key = ""
    lines = ["".join(l.split()) for l in open(f).readlines() if (l[0] != ';')]
    for l in lines:
        if l == '': continue
        if l[0] == '>':
            if curr_key != "": d[curr_key] = curr_val
            curr_key = l[1:]
            curr_val = ""
        else:
            curr_val = curr_val + l
    d[curr_key] = curr_val   
    return d

if __name__ == '__main__':

    d = fasta(sys.argv[1])
    print (d)
    


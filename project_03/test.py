from fasta import fasta
import sys

seqs = fasta(sys.argv[1])

count = {}
for key in seqs:
    for n in seqs[key]:
        count.setdefault(n, 0)
        count[n] += 1
print(count)

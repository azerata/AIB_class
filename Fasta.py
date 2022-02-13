from biopython import SeqIO

with open("seq1.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id)
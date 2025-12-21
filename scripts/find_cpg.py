import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    sys.exit("Usage: find_cpg.py input.fasta")

fasta_file = sys.argv[1]

for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    start = 0
    while True:
        idx = seq.find("CG", start)
        if idx == -1:
            break
        print(f"{record.id}\t{idx}\t{idx+2}")
        start = idx + 1

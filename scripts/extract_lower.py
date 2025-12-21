import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    sys.exit("Usage: extract_lower.py input.fasta")

fasta_file = sys.argv[1]

for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq)
    in_mask = False
    start = 0
    for i, base in enumerate(seq):
        is_lower = base.islower()
        if is_lower and not in_mask:
            start = i
            in_mask = True
        elif not is_lower and in_mask:
            print(f"{record.id}\t{start}\t{i}")
            in_mask = False
    if in_mask:
        print(f"{record.id}\t{start}\t{len(seq)}")

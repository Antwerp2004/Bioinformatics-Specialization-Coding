import numpy as np
import math
import time
from functions import *
from Bio import SeqIO
from Bio.Seq import Seq

start_time = time.perf_counter()

with open('rosalind.txt', 'r') as file:
    lines = file.readlines()
    nucleotide_1 = lines[0].strip()
    nucleotide_2 = lines[1].strip()

# records = list(SeqIO.parse('rosalind.txt', 'fasta'))
# str1 = str(records[0].seq)
# str2 = str(records[1].seq)

max_score = MiddleEdge('BLOSUM62.txt', -5, nucleotide_1, nucleotide_2, 0, len(nucleotide_1), 0, len(nucleotide_2))[1]
path = LinearSpaceAlignment('BLOSUM62.txt', -5, nucleotide_1, nucleotide_2, 0, len(nucleotide_1), 0, len(nucleotide_2))

v, w = nucleotide_1, nucleotide_2
v_idx, w_idx = 0, 0
v_aligned, w_aligned = [], []
for edge in path:
    if edge == 'D':
        v_aligned.append(v[v_idx])
        w_aligned.append(w[w_idx])
        v_idx += 1
        w_idx += 1
    elif edge == 'R':
        v_aligned.append('-')
        w_aligned.append(w[w_idx])
        w_idx += 1
    elif edge == 'Dwn':
        v_aligned.append(v[v_idx])
        w_aligned.append('-')
        v_idx += 1

with open('output.txt', 'w') as file:
    file.write(f"{max_score}\n")
    file.write(f"{"".join(v_aligned)}\n{"".join(w_aligned)}")

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"The function took {elapsed_time:.4f} seconds to run.")
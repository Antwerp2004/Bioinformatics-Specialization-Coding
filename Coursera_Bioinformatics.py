import numpy as np
import math
import time
from functions import *
from Bio import SeqIO
from Bio.Seq import Seq

start_time = time.perf_counter()

# with open('rosalind.txt', 'r') as file:
#     lines = file.readlines()
#     nucleotide_1, nucleotide_2 = lines[0].strip(), lines[1].strip()
# max_align_score, align_nucl_1, align_nucl_2 = LocalAlignment('PAM250.txt', -5, nucleotide_1, nucleotide_2)
# with open('output.txt', 'w') as file:
#     file.write(f"{max_align_score}\n")
#     file.write(f"{align_nucl_1}\n")
#     file.write(f"{align_nucl_2}")

# records = list(SeqIO.parse('rosalind.txt', 'fasta'))
# str1 = str(records[0].seq)
# str2 = str(records[1].seq)

# with open('rosalind.txt', 'r') as file:
#     lines = [line.strip() for line in file.readlines()]
#     line_0 = lines[0][1:-1].split(')(')
#     genome_P = [[int(ele) for ele in chromo.split(' ')] for chromo in line_0]
#     line_1 = lines[1][1:-1].split(')(')
#     genome_Q = [[int(ele) for ele in chromo.split(' ')] for chromo in line_1]
# two_break_dist = TwoBreakDistance(genome_P, genome_Q)
# with open('output.txt', 'w') as file:
#     file.write(f"{two_break_dist}")

with open('rosalind.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    line_0 = lines[0][1:-1].split(')(')
    genome = [[int(ele) for ele in chromo.split(' ')] for chromo in line_0]
    i1, i2, i3, i4 = [int(ele) for ele in lines[1].split(', ')]
# print(genome)
# print(f"{i1} {i2} {i3} {i4}")
two_break_genome = TwoBreakOnGenome(genome, i1, i2, i3, i4)
with open('output.txt', 'w') as file:
    for chromo in two_break_genome:
        file.write(f"({chromo[0]}") if chromo[0] < 0 else file.write(f"(+{chromo[0]}")
        for i in range(1, len(chromo)):
            file.write(f" {chromo[i]}") if chromo[i] < 0 else file.write(f"(+{chromo[i]}")
        file.write(f")")

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"The function took {elapsed_time:.4f} seconds to run.")
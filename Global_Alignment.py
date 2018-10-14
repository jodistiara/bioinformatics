# Inisialisasi Score dan Kode DNA/Asam Amino

import numpy as np
import pandas as pd

cols = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',         'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
blosum = pd.read_csv('blosum62.txt', delim_whitespace=True, header=0, index_col=0)

aa1 = "ATGC"
aa2 = "TGC"
gap = blosum['*']['*']

code1 = list("*") + list(aa1.upper())
code2 = list("*") + list(aa2.upper())

# 1. Membuat Scoring Matrix

scores = [[0 for i in range(len(code2))] for j in range(len(code1))]

# mengisi kolom pertama dan baris pertama
for i in range(len(code1)):
    scores[i][0] = gap * i
    
for j in range(len(code2)):
    scores[0][j] = gap * j

# mengisi sisanya
for i in range(1, len(code1)):
    for j in range(len(code2)):
        match = scores[i-1][j-1] + blosum[code1[i]][code2[j]]
        delete = scores[i-1][j] + gap
        insert = scores[i][j-1] + gap

        scores[i][j] = max(match, delete, insert)
        
print("Scoring Matrix: \n")
for score in scores:
    print(score)


# 2. Traceback nilai tertinggi

max_value = max(map(max, scores))
max_index = [x for x in scores if max_value in x][0]
max_index = [scores.index(max_index),max_index.index(max_value)]
print("\nMAX VALUE \nvalue: ", max_value, "\nindex: ", max_index)


# 3. Alignment Asam Amino

AlignmentA = ""
AlignmentB = ""

i = len(code1)-1
j = len(code2)-1

while (i > 0 or j > 0):
    if (i > 0) and (j > 0) and (scores[i][j] == scores[i-1][j-1] + blosum[code1[i]][code2[j]]):
        AlignmentA = code1[i] + AlignmentA
        AlignmentB = code2[j] + AlignmentB
        i = i - 1
        j = j - 1
    elif (i > 0) and (scores[i][j] == scores[i-1][j] + gap):
        AlignmentA = code1[i] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i = i - 1
    elif (j > 0) and (scores[i][j] == scores[i][j-1] + gap):
        AlignmentA = "-" + AlignmentA
        AlignmentB = code2[j] + AlignmentB
        j = j - 1

print("Amino Acid 1: ", AlignmentA)
print("Amino Acid 2: ", AlignmentB)
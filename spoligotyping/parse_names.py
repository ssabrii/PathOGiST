import re
import numpy as np
import pandas as pd

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

names = []
spoligotypes = []

with open("spo.out") as Spofile:
    for line in Spofile:
        values = line.split("\t")
        match = re.search('ERR(\d\d\d\d\d\d\d)', values[0])
        names.append(match.group(0))
        spoligotypes.append(values[1])

counter = len(names)

spoligo_distance_matrix = np.empty([counter,counter])
print(spoligo_distance_matrix)

for i in range(counter):
    for j in range(counter):
        spoligo_distance_matrix[i,j] = hamming2(spoligotypes[i], spoligotypes[j])

print(spoligo_distance_matrix)

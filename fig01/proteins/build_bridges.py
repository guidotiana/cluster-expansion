#-----------------------------------------------------------------------------#
# File: build_bridges.py
#
# Author: Matteo Tajana
# Date: 27 Sep 2022
# Description: Create bridge file starting from protein contact map and primary
# structure.
#-----------------------------------------------------------------------------#

# Import Necessary Libraries
import pdb
import glob
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Input protein pdb code
pdbCode = input("Which protein? [pdb code]\n> ")

# Create empty class protein
protein = pdb.Pdb()
print("protein:", pdbCode)
protein.ReadFile(pdbCode+".pdb", verbose=True)
protein.PrintSummary()

# Generate contact map
c = 3.7 # cutoff (in Angstrom)
print("Computing contact map for cutoff c =", c, "A...")
A = protein.GetContactMap(cutoff=c) # all atoms

# Print primary sequence (for protein.seq)
seq = protein.GetPrimaryStructure()
for l in seq:
    print(l)

# Print contacts (for protein.br)
for i in range(len(A)):
    for j in range(i,len(A)):
        if A[i][j] != 0:
            print(i+1, j+1)

# Compute graph partition function summing all Boltmann factors

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np

import statistical_potential as sp # energetics matrix

def SymLn(x):
    return np.log(np.abs(x))

# Graph parameters
M = len(sp.E) # number of configurations per node (spoiler: M=20 amino acids)
kB = 1
temperatures = np.arange(24,27,1)

grc1 = []
grc2 = []
grc3 = []
grc4 = []
grc5 = []

for T in temperatures:
    print(T)
    beta = 1./(kB*T)

    edge = 0
    path2 = 0
    path3 = 0
    clique3 = 0
    star3 = 0

    E = sp.Ecc.tolist()
    for i in range(M):
        for j in range(M):
            edge += (math.exp(-beta*E[i][j]))
            for k in range(M):
                path2 += (math.exp(-beta*(E[i][j]+E[j][k])))
                clique3 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][i])))
                for l in range(M):
                    path3 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l])))
                    star3 += (math.exp(-beta*(E[i][j]+E[i][k]+E[i][l])))


    grc1.append(edge)
    grc2.append(path2)
    grc3.append(path3)
    grc4.append(clique3)
    grc5.append(star3)
# save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4, grc5])
data = data.T
np.savetxt('1-2-3-mw.csv', data)

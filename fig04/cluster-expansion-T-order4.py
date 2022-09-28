# Compute graph contribution depending on the temperature

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np

import statistical_potential as sp # energetics matrix

# Graph parameters
M = len(sp.E) # Number of configurations per node (spoiler: M=20 amino acids)
kB = 1
temperatures = np.arange(1,300,1)

grc1 = []
grc2 = []
grc3 = []
grc4 = []
grc5 = []

for T in temperatures:
    beta = 1./(kB*T)

    edge = 0
    path2 = 0
    path3 = 0
    clique3 = 0
    star3 = 0
    paths3 = 0
    cycle4 = 0
    star4 = 0
    path4 = 0
    kite4 = 0
    Y4 = 0
    clique4 = 0
    cycle5 = 0
    star5 = 0

    E = sp.Ecc.tolist()
    for i in range(M):
        for j in range(M):
            for k in range(M):
                for l in range(M):
                    N=4
                    cycle4 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)
                    kite4 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][j])-1)
                    for m in range(M):
                        N=5
                        star4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)
                        path4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)
                        Y4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)

    grc1.append(path4)
    grc2.append(cycle4)
    grc3.append(star4)
    grc4.append(Y4)
    grc5.append(kite4)

# Save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4, grc5])
data = data.T
np.savetxt('data/4-mw.csv', data)

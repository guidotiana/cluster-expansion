# Compute graph contribution depending on the temperature

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np

# execution time measurement
from datetime import datetime
start_time = datetime.now()

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

for T in temperatures:
    # update beta
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
    path5 = 0
    cycle5 = 0
    star5 = 0
    kite5 = 0
    kite6 = 0
    cycle6 = 0
    clique5 = 0

    E = sp.Ecc.tolist()
    for i in range(M):
        for j in range(M):
            for k in range(M):
                for l in range(M):
                    N=4
                    clique4 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[k][l])-1)
                    for m in range(M):
                        N=5
                        cycle5 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][i])-1)
                        for p in range(M):
                            N=6
                            path5 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)
                            star5 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)

    grc1.append(clique4)
    grc2.append(cycle5)
    grc3.append(star5)
    grc4.append(path5)

# Save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4])
data = data.T
np.savetxt('data/5-mw.dat', data)

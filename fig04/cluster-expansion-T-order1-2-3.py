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
temperatures = np.arange(10,300,1)


grc1 = []
grc2 = []
grc3 = []
grc4 = []
grc5 = []

for T in temperatures:
    print(T)
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
            N=2
            edge += (M**(N-2))*(math.exp(-beta*E[i][j])-1)
            for k in range(M):
                N=3
                path2 += (M**(N-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)
                clique3 += (M**(N-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][i])-1)
                for l in range(M):
                    N=4
                    path3 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)
                    star3 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)

    Z1_ind = edge
    Z2_path = path2
    Z3_path = path3
    Z3_triangle = clique3
    Z3_star = star3

    grc1.append(edge)
    grc2.append(path2)
    grc3.append(path3)
    grc4.append(clique3)
    grc5.append(star3)

# Save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4, grc5])
data = data.T
np.savetxt('data/1-2-3-mw.dat', data)

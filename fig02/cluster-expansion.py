# Compute the statistical weight to the partition function of different subgraphs
# in the canonical ensemble. Slow script

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np

import statistical_potential as sp # energetics matrix

# Graph parameters
M = len(sp.E) # number of configurations per node (spoiler: M=20 amino acids)
kB = 1
T = 25
beta = 1./(kB*T)

# Graph parameters
edge = 0
path2 = 0
path3 = 0
clique3 = 0
star3 = 0
path4 = 0
cycle4 = 0
star4 = 0
Y4 = 0
kite4 = 0
clique4 = 0
path5 = 0
cycle5 = 0
star5 = 0
kite5 = 0
kite6 = 0
cycle6 = 0
path6 = 0
star6 = 0
cycle7 = 0
clique5 = 0

# Ecc = 0
E = sp.E.tolist()

for i in range(M):
    for j in range(M):
        edge += (M**(-2))*(math.exp(-beta*E[i][j])-1)
        for k in range(M):
            path2 += (M**(-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)
            clique3 += (M**(-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][i])-1)
            for l in range(M):
                path3 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)
                star3 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)
                cycle4 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)
                kite4 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][j])-1)
                clique4 += (M**(-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[k][l])-1)
                for m in range(M):
                    star4 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)
                    path4 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)
                    Y4 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)
                    cycle5 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][i])-1)
                    kite5 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)*(math.exp(-beta*E[i][m])-1)
                    clique5 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[j][m])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)*(math.exp(-beta*E[l][m])-1)
                    for p in range(M):
                        path5 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)
                        star5 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)
                        kite6 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[m][p])-1)
                        cycle6 += (M**(-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)*(math.exp(-beta*E[p][i])-1)
                        for q in range(M):
                            path6 += (M**(-7))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)*(math.exp(-beta*E[p][q])-1)
                            star6 += (M**(-7))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)*(math.exp(-beta*E[i][q])-1)


Z1_ind = edge
Z2_path = path2
Z3_path = path3
Z3_triangle = clique3
Z3_star = star3
Z4_path = path4
Z4_cycle = cycle4
Z4_star = star4
Z4_Y = Y4
Z4_kite = kite4
Z4_clique = clique4
Z5_path = path5
Z5_cycle = cycle5
Z5_star = star5
Z6_cycle = cycle6
Z6_star = star6
Z6_path = path6
Z10_clique = clique5

print('\t1-path    =', Z1_ind)
print('\t2-path    =', Z2_path)
print('\t3-path    =', Z3_path)
print('\t3-cycle   =', Z3_triangle)
print('\t3-star    =', Z3_star)
print('\t4-path    =', Z4_path)
print('\t4-cycle   =', Z4_cycle)
print('\t4-star    =', Z4_star)
print('\t4-Y       =', Z4_Y)
print('\t4-kite    =', Z4_kite)
print('\t4-clique  =', Z4_clique)
print('\t5-path    =', Z5_path)
print('\t5-cycle   =', Z5_cycle)
print('\t5-star    =', Z5_star)
print('\t6-cycle   =', Z6_cycle)
print('\t6-path    =', Z6_path)
print('\t6-star    =', Z6_star)
print('\t5-kite    =', kite5)
print('\t6-kite    =', kite6)
print('\t5-clique  =', Z10_clique)

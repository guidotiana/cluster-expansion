# Compute the statistical weight to the partition function of different subgraphs
# in the canonical ensemble. Slow script

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np
import scipy
import statistics

import statistical_potential as sp # energetics matrix
import shuffler as shf             # shuffle function

# Graph parameters
M = len(sp.E) # Number of configurations per node (spoiler: M=20 amino acids)
# Natural Units
kB = 1
T = 25
beta = 1./(kB*T)

# Perform different simulations and average the results
nsimulations = 100
print('Number of simulations =', nsimulations)
Z1_ind = []
Z2_indep = []
Z2_path = []
Z3_indep = []
Z3_path = []
Z3_triangle = []
Z3_star = []
Z4_indep = []
Z4_path_independents = []
Z4_cycle = []
Z4_kite = []
Z4_star = []
Z4_path = []
Z4_Y = []
Z5_cycle = []
Z5_path = []
Z5_star = []
Z6_clique = []

# Simulations
for n in range(nsimulations):

    # Shuffle the energy matrix...
    Eshuffle = shf.Shuffle(sp.E) # shuffled energetics matrix

    # order 1
    edge = 0
    # order 2
    path2 = 0
    # order 3
    path3 = 0
    clique3 = 0
    star3 = 0
    # order 4
    sigmapi4 = 0
    cycle4 = 0
    star4 = 0
    path4 = 0
    kite4 = 0
    Y4 = 0
    clique4 = 0
    path5 = 0
    cycle5 = 0
    star5 = 0

    # loop over the amino acids
    Es = Eshuffle.tolist()
    for i in range(M):
        for j in range(M):
            edge += math.exp(-beta*Es[i][j])-1
            for k in range(M):
                path2 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)
                clique3 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][i])-1)
                for l in range(M):
                    path3 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)
                    star3 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[i][k])-1)*(math.exp(-beta*Es[i][l])-1)
                    cycle4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[l][i])-1)
                    kite4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[l][j])-1)
                    clique4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[i][k])-1)*(math.exp(-beta*Es[i][l])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[j][l])-1)*(math.exp(-beta*Es[k][l])-1)
                    for m in range(M):
                        star4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[i][k])-1)*(math.exp(-beta*Es[i][l])-1)*(math.exp(-beta*Es[i][m])-1)
                        path4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[l][m])-1)
                        Y4 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[k][m])-1)
                        cycle5 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[l][m])-1)*(math.exp(-beta*Es[m][i])-1)
                        for p in range(M):
                            path5 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[j][k])-1)*(math.exp(-beta*Es[k][l])-1)*(math.exp(-beta*Es[l][m])-1)*(math.exp(-beta*Es[m][p])-1)
                            star5 += (math.exp(-beta*Es[i][j])-1)*(math.exp(-beta*Es[i][k])-1)*(math.exp(-beta*Es[i][l])-1)*(math.exp(-beta*Es[i][m])-1)*(math.exp(-beta*Es[i][p])-1)


    # order 1
    Z1_ind.append(edge)
    # order 2
    Z2_path.append(path2)
    # order 3
    Z3_path.append(path3)
    Z3_triangle.append(clique3)
    Z3_star.append(star3)
    # order 4
    Z4_path.append(path4)
    Z4_cycle.append(cycle4)
    Z4_kite.append(kite4)
    Z4_star.append(star4)
    Z4_Y.append(Y4)
    # order 5
    Z5_path.append(path5)
    Z5_cycle.append(cycle5)
    Z5_star.append(star5)
    # order 6
    Z6_clique.append(clique4)


print('T =', T, 'K')
print('(1,',statistics.mean(Z1_ind),      ') += (0,', statistics.stdev(Z1_ind),      ') -= (0,', statistics.stdev(Z1_ind),')')
print('(2,',statistics.mean(Z2_path),     ') += (0,', statistics.stdev(Z2_path),     ') -= (0,', statistics.stdev(Z2_path),')')
print('(3,',statistics.mean(Z3_path),     ') += (0,', statistics.stdev(Z3_path),     ') -= (0,', statistics.stdev(Z3_path),')')
print('(4,',statistics.mean(Z3_triangle), ') += (0,', statistics.stdev(Z3_triangle), ') -= (0,', statistics.stdev(Z3_triangle), ')')
print('(5,',statistics.mean(Z3_star),     ') += (0,', statistics.stdev(Z3_star),     ') -= (0,', statistics.stdev(Z3_star), ')')
print('(6,',statistics.mean(Z4_path),     ') += (0,', statistics.stdev(Z4_path),     ') -= (0,', statistics.stdev(Z4_path), ')')
print('(7,',statistics.mean(Z4_cycle),    ') += (0,', statistics.stdev(Z4_cycle),    ') -= (0,', statistics.stdev(Z4_cycle), ')')
print('(8,',statistics.mean(Z4_kite),     ') += (0,', statistics.stdev(Z4_kite),     ') -= (0,', statistics.stdev(Z4_kite), ')')
print('(9,',statistics.mean(Z4_star),     ') += (0,', statistics.stdev(Z4_star),     ') -= (0,', statistics.stdev(Z4_star), ')')
print('(10,',statistics.mean(Z4_Y),       ') += (0,', statistics.stdev(Z4_Y),        ') -= (0,', statistics.stdev(Z4_Y), ')')
print('(11,',statistics.mean(Z6_clique),  ') += (0,', statistics.stdev(Z6_clique),   ') -= (0,', statistics.stdev(Z6_clique), ')')
print('(12,',statistics.mean(Z5_path),    ') += (0,', statistics.stdev(Z5_path),     ') -= (0,', statistics.stdev(Z5_path), ')')
print('(13,',statistics.mean(Z5_cycle),   ') += (0,', statistics.stdev(Z5_cycle),    ') -= (0,', statistics.stdev(Z5_cycle), ')')
print('(14,',statistics.mean(Z5_star),    ') += (0,', statistics.stdev(Z5_star),     ') -= (0,', statistics.stdev(Z5_star), ')')

# Compute the statistical weight to the partition function of different subgraphs
# in the canonical ensemble for random gaussian energies matrix. Slow script

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np
import scipy
import statistics

import statistical_potential as sp # energetics matrix

# Perform different simulations and average the results
nsimulations = 100
print('Number of simulations =', nsimulations)

# Some parameters
# Natural units
kB = 1
T = 25
beta = 1./(kB*T)

# Data containers
Z1_ind = []
Z2_path = []
Z3_path = []
Z3_triangle = []
Z3_star = []
Z4_path = []
Z4_cycle = []
Z4_star = []
Z4_Y = []
Z4_kite = []
Z5_cycle = []
Z5_path = []
Z5_star = []
Z5_kite = []
Z6_clique = []
Z6_cycle = []
Z6_kite = []
Z10_clique = []
Zclique5 = []
Z6_star = []
Z7_cycle = []
Z7_star = []


# Statistical potentials distribution
print('Ecc matrix')
mu = sp.Ecc.mean()
std = sp.Ecc.std()*math.sqrt(2) # sqrt needed for truly corresponding symm matrix variance

# Energy extrema, uncomment if needed
#Emin = 0
#Emax = 0

counter = 0
for n in range(nsimulations):
    counter += 1
    print('simulation n', counter)
    # Random gaussian energies matrix
    Erand = np.random.normal(mu, std, (20, 20))
    # Symmetrize...
    Erand = (Erand + Erand.T)/2

    # Uncomment to fix energy extrema
    #ddd = np.diagonal(Erand)
    #if np.min(Erand) < Emin and np.min(Erand) != np.min(ddd):
    #    Emin = np.min(Erand)
    #if np.max(Erand) > Emax:
    #    Emax = np.max(Erand)
    #if np.min(ddd) < dmin:
    #    dmin = np.min(ddd)
    #print('Emin =', Emin, 'dmin =', dmin)

    M = len(Erand) # number of configurations per node (spoiler: M=20 amino acids)

    # Uncomment to fix diagonal elements to <E>
    #meanE = np.mean(Erand)
    #for i in range(M):
    #    Erand[i][i] = meanE

    # graph types
    edge = 0
    path2 = 0
    path3 = 0
    clique3 = 0
    star3 = 0
    paths3 = 0
    paths4 = 0
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
    cycle6 = 0
    star6 = 0
    cycle7 = 0

    # Loops over different amino acids
    E = Erand.tolist()

    for i in range(M):
        for j in range(M):
            edge += M**(-2)*(math.exp(-beta*E[i][j])-1)
            for k in range(M):
                path2 += M**(-3)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)
                clique3 += M**(-3)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][i])-1)
                for l in range(M):
                    path3 += M**(-4)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)
                    star3 += M**(-4)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)
                    cycle4 += M**(-4)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)
                    kite4 += M**(-4)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][j])-1)
                    clique4 += M**(-4)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[k][l])-1)
                    for m in range(M):
                        star4 += M**(-5)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)
                        path4 += M**(-5)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)
                        Y4 += M**(-5)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)
                        cycle5 += M**(-5)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][i])-1)
                        kite5 += M**(-5)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)*(math.exp(-beta*E[i][m])-1)
                        clique5 += (M**(-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[j][m])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)*(math.exp(-beta*E[l][m])-1)
                        for p in range(M):
                            path5 += M**(-6)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)
                            star5 += M**(-6)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)
                            kite6 += M**(-6)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[m][p])-1)
                            cycle6 += M**(-6)*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)*(math.exp(-beta*E[p][i])-1)
                            for q in range(M):
                                cycle7 += (M**(-7))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][p])-1)*(math.exp(-beta*E[p][q])-1)*(math.exp(-beta*E[q][i])-1)
                                star6 += (M**(-7))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)*(math.exp(-beta*E[i][q])-1)


    Z1_ind.append(edge)
    Z2_path.append(path2)
    Z3_path.append(path3)
    Z3_triangle.append(clique3)
    Z3_star.append(star3)
    Z4_cycle.append(cycle4)
    Z4_kite.append(kite4)
    Z4_star.append(star4)
    Z4_path.append(path4)
    Z4_Y.append(Y4)
    Z5_cycle.append(cycle5)
    Z5_path.append(path5)
    Z5_star.append(star5)
    Z6_clique.append(clique4)
    Z6_cycle.append(cycle6)
    Z6_star.append(star6)
    Z7_cycle.append(cycle7)
    Z5_kite.append(kite5)
    Z6_kite.append(kite6)
    Z10_clique.append(clique5)

print('T =', T, 'K')
print('(1,',statistics.mean(Z1_ind),      ') += (0,', statistics.stdev(Z1_ind),      ')')
print('(2,',statistics.mean(Z2_path),     ') += (0,', statistics.stdev(Z2_path),     ')')
print('(3,',statistics.mean(Z3_path),     ') += (0,', statistics.stdev(Z3_path),     ')')
print('(4,',statistics.mean(Z3_triangle), ') += (0,', statistics.stdev(Z3_triangle), ')')
print('(5,',statistics.mean(Z3_star),     ') += (0,', statistics.stdev(Z3_star),     ')')
print('(6,',statistics.mean(Z4_path),     ') += (0,', statistics.stdev(Z4_path),     ')')
print('(7,',statistics.mean(Z4_cycle),    ') += (0,', statistics.stdev(Z4_cycle),    ')')
print('(8,',statistics.mean(Z4_kite),     ') += (0,', statistics.stdev(Z4_kite),     ')')
print('(9,',statistics.mean(Z4_star),     ') += (0,', statistics.stdev(Z4_star),     ')')
print('(10,',statistics.mean(Z4_Y),       ') += (0,', statistics.stdev(Z4_Y),        ')')
print('(11,',statistics.mean(Z6_clique),  ') += (0,', statistics.stdev(Z6_clique),   ')')
print('(12,',statistics.mean(Z5_path),    ') += (0,', statistics.stdev(Z5_path),     ')')
print('(13,',statistics.mean(Z5_cycle),   ') += (0,', statistics.stdev(Z5_cycle),    ')')
print('(14,',statistics.mean(Z5_star),    ') += (0,', statistics.stdev(Z5_star),     ')')
print('(15,',statistics.mean(Z6_cycle),   ') += (0,', statistics.stdev(Z6_cycle),    ')')
print('(16,',statistics.mean(Z5_kite),    ') += (0,', statistics.stdev(Z5_kite),     ')')
print('(17,',statistics.mean(Z6_kite),    ') += (0,', statistics.stdev(Z6_kite),     ')')
print('(18,',statistics.mean(Z10_clique), ') += (0,', statistics.stdev(Z10_clique),  ')')
print('(19,',statistics.mean(Z6_star),    ') += (0,', statistics.stdev(Z6_star),     ')')
print('(20,',statistics.mean(Z7_cycle),   ') += (0,', statistics.stdev(Z7_cycle),    ')')

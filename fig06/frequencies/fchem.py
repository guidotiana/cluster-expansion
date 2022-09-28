# Compute the average particle number of each amino acid
# in the grandcanonical ensemble.
# Some chuncks of codes are commented so that only specific chemical potentials
# are computed.
# Pick the ones you like most.

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np
import scipy

import statistical_potential as sp # energetics matrix
import chemical_potential as cp    # amino acids chemical potentials


# Simulation parameters
M = len(sp.E) # Number of configurations per node (spoiler: M=20 amino acids)
# Natural units
kB = 1
T = 25
beta = 1./(kB*T)


# graph contributions
edge = 0
path2 = 0
clique3 = 0
star3 = 0
path3 = 0
cycle4 = 0
kite4 = 0
path4 = 0
star4 = 0
clique4 = 0
Y4 = 0
path5 = 0
cycle5 = 0
star5 = 0
cycle6 = 0
Gamma = 0 # triangle x edge

# Accessing lists is faster than numpy arrays/matrices!
E = sp.Ecc.tolist()

# mu
for i in range(M):
    for j in range(M):
        edge += (M**(-2))*(math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j])))
        for k in range(M):
            path2 += (M**(-3))*(math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k])))
            clique3 += (M**(-3))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][i]-cp.mu[i]-cp.mu[j]-cp.mu[k])))
            for l in range(M):
                star3 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l])))
                cycle4 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l])))
        '''
                path3 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          - 2*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                          - math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          + 3*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                          - 1)
                kite4 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][j]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          - math.exp(-beta*(E[i][j]+E[j][k]+E[k][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                          - math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          + 5*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                          + math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                          - 4*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                          + 1)
                clique4 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[j][k]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 2*math.exp(-beta*(E[i][j]+E[i][l]+E[j][k]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 4*math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            + 8*math.exp(-beta*(E[i][j]+E[i][l]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            + 4*math.exp(-beta*(E[i][k]+E[j][k]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            + 2*math.exp(-beta*(E[i][j]+E[i][k]+E[j][l]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            + math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 12*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 4*math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 4*math.exp(-beta*(E[i][j]+E[j][k]+E[k][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                            + 12*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                            + 3*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                            - 6*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                            + 1)
                for m in range(M):
                    path4 += (M**(-5))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                              - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                              - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                              + 3*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                              + 3*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                              - 4*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                              + 1)
                    star4 += (M**(-5))*(math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[i][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                              - 4*math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                              + 6*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                              - 4*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                              + 1)
                    Y4 += (M**(-5))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[k][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                           - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                           - math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                           - math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                           + 4*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                           + 2*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                           - 4*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                           + 1)
                    cycle5 += (M**(-5))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                               - 5*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                               + 5*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                               + 5*math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                               - 5*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                               - 5*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                               + 5*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                               - 1)
                    for p in range(M):
                        path5 += (M**(-6))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                                  - 2*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  + 3*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                                  + 6*math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                                  + math.exp(-beta*(E[i][j]+E[k][l]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - 4*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                                  - 6*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                                  + 5*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                                  - 1)
                        star5 += (M**(-6))*(math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[i][m]+E[i][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - 5*math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[i][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                                  + 10*math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                                  - 10*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                                  + 5*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                                  - 1)
                        cycle6 += (M**(-6))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][p]+E[p][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - 6*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  + 6*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                                  + 6*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  + 3*math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  - 6*math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                                  - 12*math.exp(-beta*(E[i][j]+E[j][k]+E[l][m]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]))
                                  - 2*math.exp(-beta*(E[i][j]+E[k][l]+E[m][p]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]-cp.mu[m]-cp.mu[p]))
                                  + 6*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                                  + 9*math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                                  - 6*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                                  + 1)
        '''

mZ1_edge = edge
mZ2_path = path2
mZ3_clique = clique3
mZ3_star = star3
mZ4_cycle = cycle4

mdm = []
nnn_edge = []
nnn_2path = []
nnn_3clique = []
nnn_3star = []
nnn_4cycle = []

# mu + dmu
for aa in range(20): # cycle through all amino acids
    cp.mu[aa] += cp.dmu[aa]

    edge = 0
    path2 = 0
    clique3 = 0
    star3 = 0
    cycle4 = 0
    path3 = 0

    for i in range(M):
        for j in range(M):
            edge += (M**(-2))*(math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j])))
            for k in range(M):
                path2 += (M**(-3))*(math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k])))
                clique3 += (M**(-3))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][i]-cp.mu[i]-cp.mu[j]-cp.mu[k])))
                for l in range(M):
                    star3 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l])))
                    cycle4 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][i]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l])))
                '''
                    path3 += (M**(-4))*(math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                              - 2*math.exp(-beta*(E[i][j]+E[j][k]-cp.mu[i]-cp.mu[j]-cp.mu[k]))
                              - math.exp(-beta*(E[i][j]+E[k][l]-cp.mu[i]-cp.mu[j]-cp.mu[k]-cp.mu[l]))
                              + 3*math.exp(-beta*(E[i][j]-cp.mu[i]-cp.mu[j]))
                              - 1)
                '''

    mdmZ1_edge = edge
    mdmZ2_path = path2
    mdmZ3_clique = clique3
    mdmZ3_star = star3
    mdmZ4_cycle = cycle4

    cp.mu[aa] -= cp.dmu[aa] # set default values back

    # frequencies check -> compute derivative
    naa_edge = (1/beta)*(np.log(mdmZ1_edge) - np.log(mZ1_edge))/(cp.dmu[aa])
    naa_2path = (1/beta)*(np.log(mdmZ2_path) - np.log(mZ2_path))/(cp.dmu[aa])
    naa_3clique = (1/beta)*(np.log(mdmZ3_clique) - np.log(mZ3_clique))/(cp.dmu[aa])
    naa_3star = (1/beta)*(np.log(mdmZ3_star) - np.log(mZ3_star))/(cp.dmu[aa])
    naa_4cycle = (1/beta)*(np.log(mdmZ4_cycle) - np.log(mZ4_cycle))/(cp.dmu[aa])

    nnn_edge.append(naa_edge)
    nnn_2path.append(naa_2path)
    nnn_3clique.append(naa_3clique)
    nnn_3star.append(naa_3star)
    nnn_4cycle.append(naa_4cycle)

counter = 0
frequencies = [0.083, 0.017, 0.053, 0.062, 0.039, 0.072, 0.022, 0.052, 0.057, 0.09, 0.024, 0.044, 0.051, 0.04, 0.057, 0.069, 0.058, 0.066, 0.013, 0.032]
for n in range(len(nnn_edge)):
    print('----',
          '\nmu['+str(counter)+'] =', cp.mu[counter],
          '\nreal =', frequencies[counter],
          '\nfrequency(1-pth) =', (nnn_edge[n]/np.sum(nnn_edge)),
          '\nfrequency(2-pth) =', (nnn_2path[n]/np.sum(nnn_2path)),
          '\nfrequency(3-clq) =', (nnn_3clique[n]/np.sum(nnn_3clique)),
          '\nfrequency(3-str) =', (nnn_3star[n]/np.sum(nnn_3star)),
          '\nfrequency(4-cyc) =', (nnn_4cycle[n]/np.sum(nnn_4cycle)),
         )
    counter += 1

print(np.sum(frequencies))
print(np.sum(nnn_edge))
print(np.sum(nnn_2path))
print(np.sum(nnn_3clique))
print(np.sum(nnn_3star))
print(np.sum(nnn_4cycle))

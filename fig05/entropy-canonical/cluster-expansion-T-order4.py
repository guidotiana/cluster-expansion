# Compute graph contribution depending on the temperature

import matplotlib.pyplot as plt
import math
import networkx as nx
import numpy as np

# execution time measurement
from datetime import datetime
start_time = datetime.now()

import statistical_potential as sp # energetics matrix


# Graph parameters
M = len(sp.E) # number of configurations per node (spoiler: M=20 amino acids)
kB = 1
temperatures = np.arange(24,27,1)

#sp.E += np.abs(np.amin(sp.E)) # avoid overflows

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
#            N=2 #del
#            edge += (M**(N-2))*(math.exp(-beta*E[i][j])-1)
            for k in range(M):
#                N=3 #del
#                path2 += (M**(N-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)
#                clique3 += (M**(N-3))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][i])-1)
                for l in range(M):
                    N=4 #del
#                    path3 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)
#                    star3 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)
                    cycle4 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][i])))
                    kite4 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][j])))
#                    clique4 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[j][l])-1)*(math.exp(-beta*E[k][l])-1)
                    for m in range(M):
                        N=5 #del
                        star4 += (math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[i][m])))
                        path4 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m])))
                        Y4 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[k][m])))
#                        cycle5 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)*(math.exp(-beta*E[m][i])-1)
#                        for p in range(M):
#                            N=6 #del
#                            star5 += (M**(N-6))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)*(math.exp(-beta*E[i][p])-1)


#    Z1_ind = edge
#    Z2_path = path2
#    Z3_path = path3
#    Z3_triangle = clique3
#    Z3_star = star3
#    Z4_path = path4
#    Z4_cycle = cycle4
#    Z4_star = star4
#    Z4_Y = Y4
#    Z4_kite = kite4
#    Z4_clique = clique4
#    Z5_cycle = cycle5
#    Z5_star = star5

    grc1.append(path4)
    grc2.append(cycle4)
    grc3.append(star4)
    grc4.append(Y4)
    grc5.append(kite4)

# save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4, grc5])
data = data.T
np.savetxt('Sca/4-mw.01.csv', data)


'''
plt.figure()
plt.yscale('symlog')
plt.plot(temperatures,grc1,label='4-path')
plt.plot(temperatures,grc2,label='4-cycle')
plt.plot(temperatures,grc3,label='4-star')
plt.plot(temperatures,grc4,label='4-Y')
plt.plot(temperatures,grc5,label='4-kite')
plt.xlabel(r'$T$')
plt.ylabel(r'Contribution to $Z$')
plt.title('4-edge graphs')
plt.legend()
plt.grid(True)
'''

end_time = datetime.now()
print('\nTotal execution time: {}'.format(end_time-start_time))

'''
plt.show()
'''

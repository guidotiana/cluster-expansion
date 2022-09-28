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
#                    cycle4 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][i])-1)
#                    kite4 += (M**(N-4))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][j])-1)
                    clique4 += (math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[j][k]+E[j][l]+E[k][l])))
                    for m in range(M):
                        N=5 #del
#                        star4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[i][k])-1)*(math.exp(-beta*E[i][l])-1)*(math.exp(-beta*E[i][m])-1)
#                        path4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[l][m])-1)
#                        Y4 += (M**(N-5))*(math.exp(-beta*E[i][j])-1)*(math.exp(-beta*E[j][k])-1)*(math.exp(-beta*E[k][l])-1)*(math.exp(-beta*E[k][m])-1)
                        cycle5 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][i])))
                        for p in range(M):
                            path5 += (math.exp(-beta*(E[i][j]+E[j][k]+E[k][l]+E[l][m]+E[m][p])))
                            star5 += (math.exp(-beta*(E[i][j]+E[i][k]+E[i][l]+E[i][m]+E[i][p])))


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

    grc1.append(clique4)
    grc2.append(cycle5)
    grc3.append(star5)
    grc4.append(path5)

# save numerical data
data = np.array([temperatures, grc1, grc2, grc3, grc4])
data = data.T
np.savetxt('Sca/5-mw.01.csv', data)
#np.savetxt('outdata/5-mw.02.dat', data)


'''
plt.figure()
plt.yscale('symlog')
plt.plot(temperatures,grc1,label='4-clique')
plt.plot(temperatures,grc2,label='5-cycle')
plt.plot(temperatures,grc3,label='5-star')
plt.plot(temperatures,grc4,label='5-path')
plt.xlabel(r'$T$')
plt.ylabel(r'Contribution to $Z$')
plt.title('5- and 6-edge graphs')
plt.legend()
plt.grid(True)
'''

end_time = datetime.now()
print('\nTotal execution time: {}'.format(end_time-start_time))

'''
plt.show()
'''

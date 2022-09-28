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
np.savetxt('Sca/123-mw.01.csv', data)


'''
plt.figure()
plt.yscale('symlog')
plt.plot(temperatures,grc1,label='edge')
plt.plot(temperatures,grc2,label='2-path')
plt.plot(temperatures,grc3,label='3-path')
plt.plot(temperatures,grc4,label='3-clique')
plt.plot(temperatures,grc5,label='3-star')
plt.xlabel(r'$T$')
plt.ylabel(r'Contribution to $Z$')
plt.title('1-, 2- and 3-edge Mayer factors')
plt.legend()
plt.grid(True)
'''

end_time = datetime.now()
print('\nTotal execution time: {}'.format(end_time-start_time))

'''
plt.show()
'''

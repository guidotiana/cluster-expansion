import glob
import math
import matplotlib.pyplot as plt
import numpy as np

# skip header of csv file to avoid nan
graphs_data = glob.glob('Zs-gc.csv')
data = np.genfromtxt(graphs_data[0], delimiter=' ', names=['T', 'edge', '2path', '3path', '3clique', '3star', '4path', '4cycle', '4star', '4Y', '4kite', '4clique', '5path', '5cycle', '5star'], skip_header=1)


num_lines = len(data['T'])
#dT = data['T'][1] - data['T'][0]

T = []
E_edge = []
E_edgeb = []
E_2path = []
E_3path = []
E_3clique = []
E_3star = []
E_4path = []
E_4cycle = []
E_4star = []
E_4Y = []
E_4kite = []
E_4clique = []
E_5path = []
E_5cycle = []
E_5star = []

# Note: the 'logs-vs-T-link-cc.csv' file contains the `white' Mayer factors.
# Thus, from them it is possible to compute the graph entropy calculating the whole partition function.
# In order to do so, one should sum all the `black' = `white'*(M**N) factors and cancel the similar terms.
#
# Moreover, 'logs-vs-T-link-cc.csv' has real values of fij. We must consider that
#   U = T*T*d(ln(Z))/dT
# and so exploit numpy to compute the graph entropy.

for i in range(num_lines):
    Ti = data['T'][i]
    T.append(Ti)
    if i==0:
        E_edge.append(Ti*Ti*np.log(data['edge'][0]))
        E_2path.append(Ti*Ti*np.log(data['2path'][0]))
        E_3path.append(Ti*Ti*np.log(data['3path'][0]))
        E_3clique.append(Ti*Ti*np.log(data['3clique'][0]))
        E_3star.append(Ti*Ti*np.log(data['3star'][0]))
        E_4path.append(Ti*Ti*np.log(data['4path'][0]))
        E_4cycle.append(Ti*Ti*np.log(data['4cycle'][0]))
        E_4star.append(Ti*Ti*np.log(data['4star'][0]))
        E_4Y.append(Ti*Ti*np.log(data['4Y'][0]))
        E_4kite.append(Ti*Ti*np.log(data['4kite'][0]))
        E_4clique.append(Ti*Ti*np.log(data['4clique'][0]))
        E_5path.append(Ti*Ti*np.log(data['5path'][0]))
        E_5cycle.append(Ti*Ti*np.log(data['5cycle'][0]))
        E_5star.append(Ti*Ti*np.log(data['5star'][0]))
    else:
        dT = data['T'][i] - data['T'][i-1]
        E_edge.append(Ti*Ti*(np.log(data['edge'][i])-np.log(data['edge'][i-1]))/dT)
        E_2path.append(Ti*Ti*(np.log(data['2path'][i])-np.log(data['2path'][i-1]))/dT)
        E_3path.append(Ti*Ti*(np.log(data['3path'][i])-np.log(data['3path'][i-1]))/dT)
        E_3clique.append(Ti*Ti*(np.log(data['3clique'][i])-np.log(data['3clique'][i-1]))/dT)
        E_3star.append(Ti*Ti*(np.log(data['3star'][i])-np.log(data['3star'][i-1]))/dT)
        E_4path.append(Ti*Ti*(np.log(data['4path'][i])-np.log(data['4path'][i-1]))/dT)
        E_4cycle.append(Ti*Ti*(np.log(data['4cycle'][i])-np.log(data['4cycle'][i-1]))/dT)
        E_4star.append(Ti*Ti*(np.log(data['4star'][i])-np.log(data['4star'][i-1]))/dT)
        E_4Y.append(Ti*Ti*(np.log(data['4Y'][i])-np.log(data['4Y'][i-1]))/dT)
        E_4kite.append(Ti*Ti*(np.log(data['4kite'][i])-np.log(data['4kite'][i-1]))/dT)
        E_4clique.append(Ti*Ti*(np.log(data['4clique'][i])-np.log(data['4clique'][i-1]))/dT)
        E_5path.append(Ti*Ti*(np.log(data['5path'][i])-np.log(data['5path'][i-1]))/dT)
        E_5cycle.append(Ti*Ti*(np.log(data['5cycle'][i])-np.log(data['5cycle'][i-1]))/dT)
        E_5star.append(Ti*Ti*(np.log(data['5star'][i])-np.log(data['5star'][i-1]))/dT)

head = "T edge 2path 3path 3clique 3star 4path 4cycle 4star 4Y 4kite 4clique 5path 5cycle 5star"
outdata = np.array([T, E_edge, E_2path, E_3path, E_3clique, E_3star, E_4path, E_4cycle, E_4star, E_4Y, E_4kite, E_4clique, E_5path, E_5cycle, E_5star])
#head = "T edge 2path 3path 3clique 3star 4path 4cycle 4star 4Y 4kite"
#outdata = np.array([T, E_edge, E_2path, E_3path, E_3clique, E_3star, E_4path, E_4cycle, E_4star, E_4Y, E_4kite])
outdata = outdata.T
np.savetxt('E-cc-gc.csv', outdata, delimiter=' ', header=head)

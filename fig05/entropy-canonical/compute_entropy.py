import glob
import math
import matplotlib.pyplot as plt
import numpy as np

# skip header of csv file to avoid nan
#graphs_data = glob.glob('S-cc.csv')
#data = np.genfromtxt(graphs_data[0], delimiter=' ', names=['T', 'edge', '2path', '3path', '3clique', '3star', '4path', '4cycle', '4star', '4Y', '4kite'], skip_header=1)
graphs_data = glob.glob('Zs.csv')
data = np.genfromtxt(graphs_data[0], delimiter=' ', names=['T', 'edge', '2path', '3path', '3clique', '3star', '4path', '4cycle', '4star', '4Y', '4kite', '4clique', '5path', '5cycle', '5star'], skip_header=1)


num_lines = len(data['T'])
#dT = data['T'][1] - data['T'][0]

T = []
S_edge = []
S_edgeb = []
S_2path = []
S_3path = []
S_3clique = []
S_3star = []
S_4path = []
S_4cycle = []
S_4star = []
S_4Y = []
S_4kite = []
S_4clique = []
S_5path = []
S_5cycle = []
S_5star = []

# Note: the 'logs-vs-T-link-cc.csv' file contains the `white' Mayer factors.
# Thus, from them it is possible to compute the graph entropy calculating the whole partition function.
# In order to do so, one should sum all the `black' = `white'*(M**N) factors and cancel the similar terms.
#
# Moreover, 'logs-vs-T-link-cc.csv' has real values of fij. We must consider that
#   S = -dF/dT = d(T ln(Z)) / dT
# and so exploit numpy to compute the graph entropy.

for i in range(num_lines):
    Ti = data['T'][i]
    T.append(Ti)
    if i==0:
        S_edge.append(Ti*np.log(data['edge'][0]))
        S_2path.append(Ti*np.log(data['2path'][0]))
        S_3path.append(Ti*np.log(data['3path'][0]))
        S_3clique.append(Ti*np.log(data['3clique'][0]))
        S_3star.append(Ti*np.log(data['3star'][0]))
        S_4path.append(Ti*np.log(data['4path'][0]))
        S_4cycle.append(Ti*np.log(data['4cycle'][0]))
        S_4star.append(Ti*np.log(data['4star'][0]))
        S_4Y.append(Ti*np.log(data['4Y'][0]))
        S_4kite.append(Ti*np.log(data['4kite'][0]))
        S_4clique.append(Ti*np.log(data['4clique'][0]))
        S_5path.append(Ti*np.log(data['5path'][0]))
        S_5cycle.append(Ti*np.log(data['5cycle'][0]))
        S_5star.append(Ti*np.log(data['5star'][0]))
    else:
        dT = data['T'][i] - data['T'][i-1]
        # S = -d(T lnZ) / dT
        S_edge.append((Ti*np.log(data['edge'][i]) - (Ti-dT)*np.log(data['edge'][i-1])) / dT)
        S_2path.append((Ti*np.log(data['2path'][i]) - (Ti-dT)*np.log(data['2path'][i-1])) / dT)
        S_3path.append((Ti*np.log(data['3path'][i]) - (Ti-dT)*np.log(data['3path'][i-1])) / dT)
        S_3clique.append((Ti*np.log(data['3clique'][i]) - (Ti-dT)*np.log(data['3clique'][i-1])) / dT)
        S_3star.append((Ti*np.log(data['3star'][i]) - (Ti-dT)*np.log(data['3star'][i-1])) / dT)
        S_4path.append((Ti*np.log(data['4path'][i]) - (Ti-dT)*np.log(data['4path'][i-1])) / dT)
        S_4cycle.append((Ti*np.log(data['4cycle'][i]) - (Ti-dT)*np.log(data['4cycle'][i-1])) / dT)
        S_4star.append((Ti*np.log(data['4star'][i]) - (Ti-dT)*np.log(data['4star'][i-1])) / dT)
        S_4Y.append((Ti*np.log(data['4Y'][i]) - (Ti-dT)*np.log(data['4Y'][i-1])) / dT)
        S_4kite.append((Ti*np.log(data['4kite'][i]) - (Ti-dT)*np.log(data['4kite'][i-1])) / dT)
        S_4clique.append((Ti*np.log(data['4clique'][i]) - (Ti-dT)*np.log(data['4clique'][i-1])) / dT)
        S_5path.append((Ti*np.log(data['5path'][i]) - (Ti-dT)*np.log(data['5path'][i-1])) / dT)
        S_5cycle.append((Ti*np.log(data['5cycle'][i]) - (Ti-dT)*np.log(data['5cycle'][i-1])) / dT)
        S_5star.append((Ti*np.log(data['5star'][i]) - (Ti-dT)*np.log(data['5star'][i-1])) / dT)

head = "T edge 2path 3path 3clique 3star 4path 4cycle 4star 4Y 4kite 4clique 5path 5cycle 5star"
outdata = np.array([T, S_edge, S_2path, S_3path, S_3clique, S_3star, S_4path, S_4cycle, S_4star, S_4Y, S_4kite, S_4clique, S_5path, S_5cycle, S_5star])
#head = "T edge 2path 3path 3clique 3star 4path 4cycle 4star 4Y 4kite"
#outdata = np.array([T, S_edge, S_2path, S_3path, S_3clique, S_3star, S_4path, S_4cycle, S_4star, S_4Y, S_4kite])
outdata = outdata.T
np.savetxt('S-cc.csv', outdata, delimiter=' ', header=head)

import numpy as np
import statistical_potential as sp

# Source - Proteins: Structures and Molecular Properties by Thomas E. Creighton, tab 1.1
# A C D E F G H I K L M N P Q R S T V W Y
frequencies = [0.083, 0.017, 0.053, 0.062, 0.039, 0.072, 0.022, 0.052, 0.057, 0.09, 0.024, 0.044, 0.051, 0.04, 0.057, 0.069, 0.058, 0.066, 0.013, 0.032]
M = len(frequencies) # number of amino acids
kB = 1
T = 25
beta = 1./(kB*T)

mu_highT = [-np.log(2*p)/(beta*(M*(2*M-1)-1)) for p in frequencies]

mu = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
dmu = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

for i in range(M):
    #mu[i] = (np.sum(sp.Ecc)-np.sum(sp.Ecc[i])-np.log(2*frequencies[i])/beta)/(M*(2*M-1)-1)
    mu[i] = 300*(np.sum(sp.Ecc)-np.sum(sp.Ecc[i])-np.log(2*frequencies[i])/beta)/(M*(2*M-1)-1)
    #mu[i] = -10*(frequencies[i]*np.log(frequencies[i]))/beta
    #print(mu[i])

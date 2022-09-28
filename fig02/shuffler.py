# shuffle upper (and lower) statistical potential matrix elements, keep diagonal fixed

# sources:
# https://stackoverflow.com/questions/65166640/how-can-i-permute-only-certain-entries-of-numpy-2d-array
# https://stackoverflow.com/questions/59840439/random-shuffle-of-array-but-keep-diagonal-fixed

import numpy as np
from copy import copy
from random import shuffle

import statistical_potential as sp

def Shuffle(A):
    x = A.copy()
    j, i = np.meshgrid(np.arange(x.shape[0]), np.arange(x.shape[0]))
    i, j = i.flatten(), j.flatten()
    up_i, up_j = i[i < j], j[i< j]

    elems = x[up_i, up_j]
    np.random.shuffle(elems)
    x[up_i, up_j] = elems
    x[up_j, up_i] = elems
    return x

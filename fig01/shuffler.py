# sep 13 2021
# shuffle upper (and lower) statistical potential matrix elements
# recipe: (1) select upper triangle exploiting rows
#         (2) shuffle them keeping zeros on left (shape...)
#         (3) join the arrays

# sources:
# https://stackoverflow.com/questions/65166640/how-can-i-permute-only-certain-entries-of-numpy-2d-array
# https://stackoverflow.com/questions/59840439/random-shuffle-of-array-but-keep-diagonal-fixed

import numpy as np
from copy import copy
from random import shuffle

E = np.array([[-13,0,12,26,3,-7,34,-22,14,-1,25,28,10,8,43,-6,-9,-10,-9,9],
[0,-106,3,69,-23,-8,-19,16,71,-8,19,13,0,5,24,-2,19,6,8,4],
[12,3,4,-15,39,-22,-39,59,-76,67,65,-30,4,-17,-72,-31,-29,58,24,0],
[26,69,-15,-3,27,25,-45,35,-97,43,44,-32,-10,-17,-74,-26,0,34,29,-10],
[3,-23,39,27,-44,-38,-16,-19,44,-30,-42,18,20,-29,41,29,31,-22,-16,0],
[-7,-8,-22,25,-38,-38,20,25,11,23,19,-14,-11,-6,-4,-16,-26,16,18,14],
[34,-19,-39,-45,-16,20,-29,49,22,16,99,-24,-21,-2,-12,-5,-19,19,-12,-34],
[-22,16,59,35,-19,25,49,-22,36,-41,-28,53,25,36,42,21,14,-25,2,11],
[14,71,-76,-97,44,11,22,36,25,19,0,-33,11,-38,75,-13,-9,44,22,-21],
[-1,-8,67,43,-30,23,16,-41,19,-27,-20,30,42,26,35,25,20,-29,-9,24],
[25,19,65,44,-42,19,99,-28,0,-20,4,8,-34,46,31,14,19,-14,-67,-13],
[28,13,-30,-32,18,-14,-24,53,-33,30,8,-53,-18,-25,-14,-14,-11,50,6,-20],
[10,0,4,-10,20,-11,-21,25,11,42,-34,-18,26,-42,-38,1,-7,9,-28,-33],
[8,5,-17,-17,-29,-6,-2,36,-38,26,46,-25,-42,29,-52,-14,-14,24,8,-20],
[43,24,-72,-74,41,-4,-12,42,75,35,31,-14,-38,-52,11,17,-35,30,-16,-25],
[-6,-2,-31,-26,29,-16,-5,21,-13,25,14,-14,1,-14,17,-20,-8,18,34,9],
[-9,19,-29,0,31,-26,-19,14,-9,20,19,-11,-7,-14,-35,-8,3,25,22,13],
[-10,6,58,34,-22,16,19,-25,44,-29,-14,50,9,24,30,18,25,-29,-7,2],
[-9,8,24,29,-16,18,-12,2,22,-9,-67,6,-28,8,-16,34,22,-7,-12,-4],
[9,4,0,-10,0,14,-34,11,-21,24,-13,-20,-33,-20,-25,9,13,2,-4,-6]])
E = E/100

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

Eshuffled = Shuffle(E)
np.savetxt('mj_float_shuffled.mat', Eshuffled, fmt='%.2f')

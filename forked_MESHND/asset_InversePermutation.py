import numpy as np

def inverse_permutation(perm):
    inverse = [0] * len(perm)
    for i, p in enumerate(perm):
        inverse[p] = i
    return inverse

p = [0,2,3,1,8,12,10,11,14,15,9,13,4,5,6,7] # Gotta start at zero in python and C
# pInv = [0,3,1,2,12,13,14,15,4,10,6,7,5,11,8,9] # Expected answer

print('p = '+str(p))
print('inv_p = '+str(inverse_permutation(p)))

# for i, a in enumerate(p):
    # print i,a
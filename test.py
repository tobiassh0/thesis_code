from func_load import *

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

n = []
size = (100,300)
n1 = np.ones(size)
n2 = n1*2
n3 = n1*3

n.append(n1)
n.append(n2)
n.append(n3)
n = np.array(n)
print(n.shape)
m0 = np.mean(n,axis=0)
m1 = np.mean(n,axis=1)
m2 = np.mean(n,axis=2)

print(m0.shape)
print(m1.shape)
print(m2.shape)
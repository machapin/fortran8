import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import glob

filepath = 'data/ibm21_newton/z_00000001.d'
# filepath = 'data/data.d'
data_d = np.loadtxt(filepath)

print(data_d)
print(data_d.shape)

filepath = 'data/ibm21_newton/z_00000001.bin'
# filepath = 'data/data.bin'

with open(filepath, mode='rb') as file:
    data = np.fromfile(file, dtype=np.float64)

data = data.reshape(-1, 18)

print(data)
print(data.shape)

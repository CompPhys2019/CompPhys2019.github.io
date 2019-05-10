#!/usr/env/bin python

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('corr.dat')

x = np.arange(data.shape[1])

for i in range(data.shape[0]):
    plt.plot(x,data[i],label=str(i))
    plt.legend()

plt.show()
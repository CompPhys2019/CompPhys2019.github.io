#!/usr/bin/env python

import numpy as np

g = lambda a,x: np.arccos(a*x)
g = lambda a,x: (np.cos(x)+x)/(a+1)
# g = lambda a,x: (np.cos(x)-x)/(a-1)

x0 = -4
# x0 = -1.85
x0 = 1.34
for i in range(1000):
    x0 = g(1./6,x0)

print(x0)
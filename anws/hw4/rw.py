#!/usr/bin/env python

import numpy as np
import random

move = np.array([[-1,0],[1,0],[0,-1],[0,1]])

def Simple_RW(r):
    i = random.randint(0,3)
    r = r + move[i] 
    return r

Nsam = 100000

mean = np.array([0,0])
for i in range(Nsam):
    r = np.array([0,0])
    for j in range(6):
        r = Simple_RW(r)

    mean += r

print(mean/Nsam)

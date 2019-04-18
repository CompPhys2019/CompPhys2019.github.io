#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def iters(n,r,x0=0.5):
    xs=[x0]
    f = lambda x: 4*r*x*(1-x)
    for i in range(1,n): xs.append(f(xs[-1]))
    return xs

def bifurcation(x0s,rs,niters=1000,last=100,alpha=0.125):
    # x0s, sequences of x0 to be choosed
    # rs, sequences of r to be choosed
    # niters, number of iteration times
    # plot last iteration times
    for x0 in x0s:
        fig = plt.figure()  
        for r in rs:
            xn = iters(niters,r,x0)[-last:]
            r1 = last*[r]
            plt.scatter(r1,xn,marker='o',s=0.1,alpha=alpha,color='black')
        
        plt.title("Bifurcation diagram with $x_0 = {}$".format(x0))
        plt.xlabel('$r$')
        plt.ylabel('iterated values of x')
        plt.savefig('bifurcation{:.8f}.eps'.format(x0))
    plt.show()

x0s = [0.01,0.1,0.5,0.9]
# x0s = [0.9]
rs = np.arange(0.7,1.0,0.0001)

bifurcation(x0s,rs,2000,200,0.05)


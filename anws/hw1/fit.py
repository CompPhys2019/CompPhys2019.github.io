#!/usr/bin/env python

import numpy as np
from lmfit import Model
from scipy.stats import chi2
import math

data = np.loadtxt('temp.dat')

t, y1, y2 = data[:,0], data[:,2], data[:,4] 
ye = 0.1*np.sqrt(2)
# ye is error of T-Ts

# exponential fitting
def expfun(t,c,r):
    return c*np.exp(-r*t)

expm = Model(expfun, independent_vars=['t'])
expm.set_param_hint('c',value=1,vary=True)
expm.set_param_hint('r',value=1,vary=True)

result = expm.fit(y1, t=t, weights=1/ye)
print(result.fit_report(show_correl=False,sort_pars=True))
goodness=chi2.sf(result.chisqr,result.nfree)


# linear fitting
def linfun(t,c,r):
    return -r*t + c

linm = Model(linfun, independent_vars=['t'])
linm.set_param_hint('c',value=1,vary=True)
linm.set_param_hint('r',value=1,vary=True)

result = expm.fit(np.log(y1), t=t, weights=y1/ye)
print(result.fit_report(show_correl=False,sort_pars=True))
goodness=chi2.sf(result.chisqr,result.nfree)


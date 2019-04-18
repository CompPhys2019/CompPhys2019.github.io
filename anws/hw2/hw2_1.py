#!/usr/bin/env python

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

x = sym.symbols('x')
fx = x*sym.exp(x)
f = sym.lambdify(x, fx, 'numpy')
dfn = lambda n: sym.lambdify(x, fx.diff(x,n), 'numpy')(1)

diff_foward1 = lambda f,x,h: (f(x+h)-f(x))/h
diff_foward2 = lambda f,x,h: (f(x+2*h) - 2*f(x+h) + f(x))/h**2

diff_central1 = lambda f,x,h: (f(x+h)-f(x-h))/h/2
diff_central2 = lambda f,x,h: (f(x+h) - 2*f(x) + f(x-h))/h**2

diff_Richardson1 = lambda f,x,h: (f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/h/12
diff_Richardson2 = lambda f,x,h: (-f(x-2*h) + 16*f(x-h) - 30*f(x) + 16*f(x+h) - f(x+2*h))/h**2/12

h = np.arange(0.05,0.55,0.05)

err1f = np.abs(diff_foward1(f,1,h)-dfn(1))
err2f = np.abs(diff_foward2(f,1,h)-dfn(2))

err1c = np.abs(diff_central1(f,1,h)-dfn(1))
err2c = np.abs(diff_central2(f,1,h)-dfn(2))

err1r = np.abs(diff_Richardson1(f,1,h)-dfn(1))
err2r = np.abs(diff_Richardson2(f,1,h)-dfn(2))

z1f = np.polyfit(np.log(h), np.log(err1f), 1)
z1c = np.polyfit(np.log(h), np.log(err1c), 1)
z1r = np.polyfit(np.log(h), np.log(err1r), 1)
z2f = np.polyfit(np.log(h), np.log(err2f), 1)
z2c = np.polyfit(np.log(h), np.log(err2c), 1)
z2r = np.polyfit(np.log(h), np.log(err2r), 1)
print(z1f[0],z1c[0],z1r[0],z2f[0],z2c[0],z2r[0])

# print(np.log(err2f))
# print(np.log(err2c))
# print(np.log(err2r))

plt.plot(h,err1f,label='first forward')
plt.plot(h,err1c,label='first central')
plt.plot(h,err1r,label='first Richardson')
plt.plot(h,err2f,label='second forward')
plt.plot(h,err2c,label='second central')
plt.plot(h,err2r,label='second Richardson')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

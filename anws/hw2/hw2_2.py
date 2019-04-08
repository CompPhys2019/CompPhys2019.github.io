#!/usr/bin/env python

import sympy as sym

x = sym.symbols('x')

fx = sym.exp(x)/(sym.sin(x)**3+sym.cos(x)**3)



df0 = lambda n: sym.lambdify(x, fx.diff(x,n), 'numpy')(0)

print(df0(1),df0(2),df0(3),df0(4),df0(5))

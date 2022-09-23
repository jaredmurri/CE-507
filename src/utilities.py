from random import vonmisesvariate
import unittest
import math
import numpy as np
import sympy
import scipy
import matplotlib.pyplot as plt


def taylorExpansion(fun,a,order):
    x=list(fun.atoms(sympy.Symbol))[0]
    t=0
    for i in range(0,order+1):
        df=sympy.diff(fun,x,i)
        term=(df.subs(x,a)/sympy.factorial(i))*(x-a)**i
        t+=term
    return t 

def plot_taylor_sin():
    x = np.linspace(-1, 1)
    for degree in np.array([0,1,3,5,7]):
        sin_taylor = scipy.interpolate.approximate_taylor_polynomial(np.sin, 0, degree,1)
        plt.plot(x, sin_taylor(np.pi*x))
    plt.plot(x,np.sin(np.pi*x),color="black")
    plt.show()

plot_taylor_sin()
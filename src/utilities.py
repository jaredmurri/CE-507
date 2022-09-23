from random import vonmisesvariate
import unittest
import math
import numpy
import sympy
import matplotlib.pyplot as plt


def taylorExpansion(fun,a,order):
    x=list(fun.atoms(sympy.Symbol))[0]
    t=0
    for i in range(0,order+1):
        df=sympy.diff(fun,x,i)
        term=(df.subs(x,a)/sympy.factorial(i))*(x-a)**i
        t+=term
    return t 

def plottest():

    x = numpy.linspace(-1, 1)
    y = numpy.sin(numpy.pi*x)   
    plt.plot(x,y)
    plt.show()

def lagrangetry():
    from scipy.interpolate import lagrange

    f = lagrange(x, y)

    fig = plt.figure(figsize = (10,8))
    plt.plot(x_new, f(x_new), 'b', x, y, 'ro')
    plt.title('Lagrange Polynomial')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
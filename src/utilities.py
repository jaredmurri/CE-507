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

    x = numpy.linspace(0, 2 * numpy.pi, 200)
    y = numpy.sin(x)

    fig, ax = plt.subplots()
    ax.plot(x, y)
    plt.show()
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

def plot_taylor_exp():
    x = np.linspace(-1, 1)
    for degree in np.array([1,2,3,4]):
        exp_taylor = scipy.interpolate.approximate_taylor_polynomial(np.exp, 0, degree,1)
        plt.plot(x, exp_taylor(x))
    plt.plot(x,np.exp(x),color="black")
    plt.show()

def plot_taylor_erfc():
    x = np.linspace(-2, 2)
    for degree in [1,3,5,7]:
        erfc_taylor = scipy.interpolate.approximate_taylor_polynomial(scipy.special.erfc, 0, degree,1)
        plt.plot(x, erfc_taylor(x))
    plt.plot(x,scipy.special.erfc(x),color="black")
    plt.show()

def plot_monomial_basis():
    x = np.linspace(0, 1)
    for degree in range(0,11):
        plt.plot(x, x**degree)
    plt.show()


def sin_error():
    # x = np.linspace(-1, 1)
    error=[]
    degree_list = np.array([0,1,3,5,7])
    for degree in degree_list:
        sin_taylor = scipy.interpolate.approximate_taylor_polynomial(np.sin, 0, degree,1)
        err_fun = lambda x : abs( sin_taylor(np.pi * x) - np.sin(np.pi*x) )  
        error.append( scipy.integrate.quad( err_fun, -1, 1 )[0] )
    plt.plot(degree_list, error, color="black")
    plt.yscale("log")
    plt.show()
   
def exp_error():
    # x = np.linspace(-1, 1)
    error=[]
    degree_list = np.array([1,2,3,4])
    for degree in degree_list:
        exp_taylor = scipy.interpolate.approximate_taylor_polynomial(np.exp, 0, degree,1)
        err_fun = lambda x : abs( exp_taylor(x) - np.exp(x) )  
        error.append( scipy.integrate.quad( err_fun, -1, 1 )[0] )
    plt.plot(degree_list, error, color="black")
    plt.yscale("log")
    plt.show()

def erfc_error():
    # x = np.linspace(-1, 1)
    error=[]
    degree_list = np.array([1,3,5,7])
    for degree in degree_list:
        exp_taylor = scipy.interpolate.approximate_taylor_polynomial(scipy.special.erfc, 0, degree,1)
        err_fun = lambda x : abs( exp_taylor(x) - scipy.special.erfc(x) )  
        error.append( scipy.integrate.quad( err_fun, -2, 2 )[0] )
    plt.plot(degree_list, error, color="black")
    plt.yscale("log")
    plt.show()
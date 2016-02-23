#!/usr/bin/env python

# This .py file is to calculate the Nth Fibonacci number.

# import numpy
from numpy import *

# create a function
def Fibonacci(num):
    Fn=((1+sqrt(5))**num-(1-sqrt(5))**num)/2**num/sqrt(5)
    print "the %d Fibonacci number is %d."%(num,Fn)
    return Fn

# call this function
num=raw_input("please enter a number:")
num=int(num)
Fibonacci(num)
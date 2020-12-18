# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:44:12 2020

@author: Hugo
"""

import matplotlib.pyplot as plt
import numpy as np

delta_t=1
delta_r=1
A=1
B=1
C=1

Tmax=1000
rmax=1000

w=1 # injection

def iter_h(h,n,i):
    """i : indice spatial, n : indice temporel"""
    
    h[n,i+1]=h[n,i]+C**2+delta_t/(n*delta_r**2)*(h[n+1,i]**3*(n+1)*(h[n+1,i]-h[n,i])-h[n,i]**3*n*(h[n,i]-h[n-1,i]))+w
    

def iter_theta(theta,n,i):
    """i : indice spatial, n : indice temporel"""

    theta[n,i+1]=(1/h[n+1,i]**3)*(theta[n,i]*h[n,i]**3+(1/B)*(A*(delta_t/(n*delta_r**2))*((n+1)*(h[n+1,i]-h[n,i])*theta[n+1,i]*h[n+1,i]**5-n*(h[n,i]-h[n-1,i])*theta[n,i]*h[n,i]**5)-C*h[n,i]*theta[n,i]))
    

# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:44:12 2020

@author: Hugo
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


delta_t=0.1
delta_r=1
A=1
B=1
C=0.2

t_max=50
r_max=50

w=0 # injection
i_inj=10


def iter_h(h,n,i):
    """i : indice spatial, n : indice temporel"""
    # test=0
    # if test:
    #     h[n+1,i]=h[n-1,i]
    # else:
    if i>=r_max-1:
        h[n+1,i]=0
    elif i >= i_inj and i<r_max-1 :
        
        h[n+1,i]=h[n-1,i]+(C**2 *delta_t)/(4*i*delta_r**2)*((i+1/2)*(h[n,i+1]+h[n,i])**3*(h[n,i+1]-h[n,i])-(i-1/2)*(h[n,i]+h[n,i-1])**3*(h[n,i]-h[n,i-1]))

    elif i!=0 and i<r_max-1:
        h[n+1,i]=h[n-1,i]+(C**2 *delta_t)/(4*i*delta_r**2)*((i+1/2)*(h[n,i+1]+h[n,i])**3*(h[n,i+1]-h[n,i])-(i-1/2)*(h[n,i]+h[n,i-1])**3*(h[n,i]-h[n,i-1]))+2*delta_t*w

    elif i==0 and i<r_max-1:
        h[n+1,i]=h[n+1,i+1]    

def iter_theta(theta,n,i):
    """i : indice spatial, n : indice temporel"""

    theta[n+1,i]=(1/h[n,i+1]**3)*(theta[n,i]*h[n,i]**3+(1/B)*(A*(delta_t/(n*delta_r**2))*((i+1)*(h[n,i+1]-h[n,i])*theta[n,i+1]*h[n,i+1]**5-n*(h[n,i]-h[n,i-1])*theta[n,i]*h[n,i]**5)-C*h[n,i]*theta[n,i]))
    

h=np.zeros((t_max,r_max))

# Mettre profil en h(r)=h_0(t)(1-r^2/R^2)^{1/3

for i in range(0,30):
    h[0,i]=1
    h[1,i]=1



for n in range(1,t_max-2):
    for i in range(r_max-1,-1,-1):
        #print(i)
        iter_h(h, n, i)
        
# plt.plot(h[2,:])
# plt.show()       
        
fig, ax = plt.subplots()
xdata, ydata = [], []
line, = ax.plot([],[], color='blue')
ax.set_xlim([0, r_max])
ax.set_ylim([0,np.max(h)*1.1])

        


def animate(k): 
    line.set_data(np.arange(r_max),h[k,:])
    return line,


ani = animation.FuncAnimation(fig=fig, func=animate, frames=np.arange(t_max), blit=True, interval=.01, repeat=True)

plt.show()


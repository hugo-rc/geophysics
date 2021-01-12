# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:44:12 2020

@author: Hugo
"""

# =============================================================================
# Modules
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.integrate import simps

# =============================================================================
# Parametres
# =============================================================================

profil_depart=False # True pour partir d'un profil initial, False pour partir sans matiere 
i_init=20 # Rayon du profil de depart, utile seulement si profil_depart=True
h0=1 # Hauteur du profil de depart, utile seulement si profil_depart=True

kappa = 10**(-6)
delta_r = 1
delta_t = delta_r**2/(2*h0**3)
rho0 = 3.0
rhom = 3.3
drho = abs(rho0 - rhom)/rho0
T_ext=462 # temperature surface venus (degree celsius)

A = 11/20 - 11/(4*drho) + 23/(2*drho**2) - 11/(2*drho**3) - 19/(4*drho**4) - 11/(20*drho**5)
B = 2*((1/3)*(1+1/drho**3)+(1/drho)*(1+1/drho))
C = 1 + 1/drho

i_max=200
n_max=10000
t_max=100
r_max=200

w_inj=0 # vitesse injection (>1e-4)
rho_inj=rho0 
theta_inj=0
i_inj=20

# =============================================================================
# Fonctions
# =============================================================================
        
def iter_h(h,n,i):
    """i : indice spatial, n : indice temporel"""
    
    if i!=0 and i >= i_inj and i<i_max-1 :
        h[n+1,i]=h[n,i]+((C**2 *delta_t)/(2*i*delta_r**2))*(((i+1)*h[n,i+1]**3+i*h[n,i]**3)*(h[n,i+1]-h[n,i])-(i*h[n,i]**3+(i-1)*h[n,i-1]**3)*(h[n,i]-h[n,i-1]))

    elif i!=0 and i<i_inj and i<i_max-1:
        h[n+1,i]=h[n,i]+(C**2 *delta_t)/(2*i*delta_r**2)*(((i+1)*h[n,i+1]**3+i*h[n,i]**3)*(h[n,i+1]-h[n,i])-(i*h[n,i]**3+(i-1)*h[n,i-1]**3)*(h[n,i]-h[n,i-1]))+2*delta_t*w_inj

    elif i==0 :
        h[n+1,i]=h[n+1,i+1]  

def iter_theta(theta,n,i):
    """i : indice spatial, n : indice temporel
    
    WIP
    """
    if i!=0 and i >= i_inj and i<i_max-1 :

        theta[n+1,i]=(1/h[n,i+1]**3)*((theta[n,i]*h[n,i]**3)+(A*delta_t/(B*n*delta_r**2))*((h[n,i+1]-h[n,i])*((i+1)*theta[n,i+1]*h[n,i+1]**5+i*theta[n,i]*h[n,i]**5)-(h[n,i]-h[n,i-1])*(i*theta[n,i]*h[n,i]**5+(i-1)*theta[n,i-1]*h[n,i-1]**5))-C*h[n,i]*theta[n,i])
    elif i!=0 and i<i_inj and i<i_max-1:
        theta[n+1,i]=(1/h[n,i+1]**3)*((theta[n,i]*h[n,i]**3)+(A*delta_t/(B*n*delta_r**2))*((h[n,i+1]-h[n,i])*((i+1)*theta[n,i+1]*h[n,i+1]**5+i*theta[n,i]*h[n,i]**5)-(h[n,i]-h[n,i-1])*(i*theta[n,i]*h[n,i]**5+(i-1)*theta[n,i-1]*h[n,i-1]**5))-C*h[n,i]*theta[n,i])+rho_inj*theta_inj*w_inj

    elif i==0 :
        theta[n+1,i]=theta[n+1,i+1]
    
def profil_init(h,i_init):
    """ Pour partir d'un profil initial, utile notament si on travail a volume constant"""
    for i in range(0,i_init):
        r=i*delta_r
        R=i_init*delta_r
        h[0,i]=h0*(1-r**2/R**2)**(1/3)
        h[1,i]=h0*(1-r**2/R**2)**(1/3)

    return h   

def hauteur_centre(h):
    """Hauteur centre du panache au fil du temps"""
    h_centre=h[:,0]
    return h_centre

def rayon_max(h):
    """ Rayon max du panache au fil du temps"""
    shape=np.shape(h)
    R_max=[]
    for n in range(shape[0]):
        for i in range(shape[1]):
            if h[n,i]==0:
                R_max.append((i-1)*delta_r)
                break
    return R_max

def volume(h,n):
    """ Volume du panache a l'instant n"""
    r=np.arange(i_max)*delta_r
    z=h[n,:]
    V=simps(r*z,r)   
    return V


# =============================================================================
# Programme principal
# =============================================================================

h=np.zeros((n_max,i_max))
print('Taille h :', np.shape(h))


if profil_depart:
    
    h=profil_init(h,i_init)


for n in range(1,n_max-2):
    for i in range(i_max-1,-1,-1):

        iter_h(h, n, i)
        
      
        
fig, ax = plt.subplots()
xdata, ydata = [], []
line, = ax.plot([],[], color='blue')
ax.set_xlim([0, i_max])
ax.set_ylim([0,np.max(h)*1.1])
time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
time_text.set_text('t=0')
        


def animate(k): 
    line.set_data(np.arange(len(h[k,:])),h[k,:])
    time_text.set_text('time = %.1d' % k)
    return line, time_text,


ani = animation.FuncAnimation(fig=fig, func=animate, frames=np.arange(n_max), blit=True, interval=.001, repeat=False)

plt.show()


h_centre=hauteur_centre(h)
r_max=rayon_max(h)
t=np.arange(n_max)*delta_t

plt.loglog(t,h_centre)
plt.title('Hauteur centre')
plt.show()
plt.loglog(t,r_max)
plt.title('Rayon max')
plt.show()


plt.plot(np.arange(n_max),[volume(h,n) for n in range(n_max)])
plt.title('Volume total')
plt.show()

theta=np.zeros((n_max,i_max))
for n in range(1,n_max-2):
    for i in range(i_max-1,-1,-1):

        iter_theta(h, n, i)

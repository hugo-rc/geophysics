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
#import matplotlib.animation as animation
from matplotlib import animation#, rc
from scipy.integrate import simps

# =============================================================================
# Parametres
# =============================================================================

profil_depart=True # True pour partir d'un profil initial, False pour partir sans matiere 
injection=True
record = True
i_init=2 # Rayon du profil de depart, utile seulement si profil_depart=True
h0=0.24 # Hauteur du profil de depart, utile seulement si profil_depart=True

kappa = 10**(-6)
Q=100
mu=1e20 #1e19 - 1e20
g=8.87


rho0 = 4e3
rhom = 3.3e3
drho = abs(rho0 - rhom)/rho0
T_ext=462 # temperature surface venus (degree celsius)

A = 11/20 - 11/(4*drho) + 23/(2*drho**2) - 11/(2*drho**3) - 19/(4*drho**4) - 11/(20*drho**5)
B = 2*((1/3)*(1+1/drho**3)+(1/drho)*(1+1/drho))
C = 1 + 1/drho


delta_r = 1e-3
delta_t = delta_r**2/(2*h0**3*C**2)

i_max=100
n_max=1000
#t_max=100
r_max=i_max*delta_r


H=((3*mu*Q)/(rho0*g))**(1/4)
tau=((3*mu*Q)/(rho0*g*kappa**2))**(1/2)
R=((3*mu*Q**5)/(kappa**4*rho0*g))**(1/8)

rho_inj=0
theta_inj=0
i_inj=0
w_inj=0

if injection:
    #w_inj=1 # vitesse injection (>1e-4)
    rho_inj=rho0 
    theta_inj=1
    i_inj=2
    w_inj=Q*tau/H*(1/(2*np.pi*(i_inj*delta_r*R)**2))


#print(w_inj)
# =============================================================================
# Fonctions
# =============================================================================
  

def iter_h(h,n,i):
    """i : indice spatial, n : indice temporel"""
    if profil_depart:
        delta_t=delta_r**2/(2*np.max(h[n,:])**3*C**2)
    else:
        delta_t=1e-3
        
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
    
    if profil_depart:
        delta_t=delta_r**2/(2*np.max(h[n,:])**3*C**2)
    else:
        delta_t=1e-3
    seuil=0.1
    if i!=0 and i >= i_inj and i<i_max-1 : 
        if h[n,i+1]<seuil: # cas ou il n'y a pas encore depanache en i, temperature = temperature surface, +on ne résout trop proche des bords sinon ça diverge
            theta[n+1,i]=0
        else:
            tmp1=1/h[n,i+1]**3
            tmp2=theta[n,i]*h[n,i]**3
            tmp3=h[n,i+1]-h[n,i]
            tmp4=(i+1)*theta[n,i+1]*h[n,i+1]**5
            tmp5=i*theta[n,i]*h[n,i]**5
            tmp6=h[n,i]-h[n,i-1]
            tmp7=(i-1)*theta[n,i-1]*h[n,i-1]**5
            tmp8=h[n,i]*theta[n,i]
            theta[n+1,i]=tmp1*(tmp2+(A*delta_t/(B*i*delta_r**2))*(tmp3*(tmp4+tmp5)-tmp6*(tmp5+tmp7))-(C*delta_t/B)*tmp8)
#            theta[n+1,i]=(1/h[n,i+1]**3)*((theta[n,i]*h[n,i]**3)+(A*delta_t/(B*i*delta_r**2))*((h[n,i+1]-h[n,i])*((i+1)*theta[n,i+1]*h[n,i+1]**5+i*theta[n,i]*h[n,i]**5)-(h[n,i]-h[n,i-1])*(i*theta[n,i]*h[n,i]**5+(i-1)*theta[n,i-1]*h[n,i-1]**5))-(C*delta_t/B)*h[n,i]*theta[n,i])
    elif i!=0 and i<i_inj and i<i_max-1 : 
        if h[n,i+1]<seuil: # cas ou il n'y a pas encore depanache en i, temperature = temperature surface
            theta[n+1,i]=0
        else:
            tmp1=1/h[n,i+1]**3
            tmp2=theta[n,i]*h[n,i]**3
            tmp3=h[n,i+1]-h[n,i]
            tmp4=(i+1)*theta[n,i+1]*h[n,i+1]**5
            tmp5=i*theta[n,i]*h[n,i]**5
            tmp6=h[n,i]-h[n,i-1]
            tmp7=(i-1)*theta[n,i-1]*h[n,i-1]**5
            tmp8=h[n,i]*theta[n,i]
            theta[n+1,i]=tmp1*(tmp2+(A*delta_t/(B*i*delta_r**2))*(tmp3*(tmp4+tmp5)-tmp6*(tmp5+tmp7))-(C*delta_t/B)*tmp8)+rho_inj*theta_inj*w_inj
#            theta[n+1,i]=(1/h[n,i+1]**3)*((theta[n,i]*h[n,i]**3)+(A*delta_t/(B*i*delta_r**2))*((h[n,i+1]-h[n,i])*((i+1)*theta[n,i+1]*h[n,i+1]**5+i*theta[n,i]*h[n,i]**5)-(h[n,i]-h[n,i-1])*(i*theta[n,i]*h[n,i]**5+(i-1)*theta[n,i-1]*h[n,i-1]**5))-(C*delta_t/B)*h[n,i]*theta[n,i])+rho_inj*theta_inj*w_inj

    elif i==0 :
        theta[n+1,i]=theta[n+1,i+1]
        
    
def profil_init(h,i_init):
    """ Pour partir d'un profil initial, utile notament si on travail a volume constant"""
    for i in range(0,i_init):
        r=i*delta_r
        R=i_init*delta_r
        h[0,i]=h0*(1-r**2/R**2)**(1/3)
#        h[1,i]=h0*(1-r**2/R**2)**(1/3)

    return h   

def hauteur_centre(h):
    """Hauteur centre du panache au fil du temps"""
    h_centre=h[:,0]
    return np.array(h_centre)

def rayon_max_instant(h,n):
    """Indice rayon max du panache a un instant n"""
    shape=np.shape(h)
    for i in range(shape[1]):
            if h[n,i]==0:
                return i-1
                

def rayon_max(h):
    """ Rayon max du panache au fil du temps"""
    shape=np.shape(h)
    R_max=[]
    for n in range(shape[0]):
        for i in range(shape[1]):
            if h[n,i]==0:
                R_max.append((i-1)*delta_r)
                break
    return np.array(R_max)

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


for n in range(0,n_max-1):
    for i in range(i_max-1,-1,-1):

        iter_h(h, n, i)
        
h_centre=hauteur_centre(h)
r_max=rayon_max(h)
t=np.arange(n_max)*(delta_r**2/(2*h_centre**3*C**2))      
        
fig, ax = plt.subplots()
xdata, ydata = [], []
line, = ax.plot([],[], color='blue')
ax.set_xlabel('r/R')
ax.set_ylabel('h/H')
ax.set_xlim([0, i_max*delta_r])
ax.set_ylim([0,np.max(h)*1.1])
time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
time_text.set_text('t/tau=0')
        


def animate(k): 
    line.set_data(np.arange(len(h[k,:]))*delta_r,h[k,:])
    time_text.set_text('frame = %.1d' % k)
    return line, time_text,


ani = animation.FuncAnimation(fig=fig, func=animate, frames=np.arange(n_max), blit=True, interval=.001, repeat=False, save_count=n_max)

if record:
    plt.rcParams['animation.ffmpeg_path'] =r'E:\Anaconda\pkgs\ffmpeg-4.3.1-ha925a31_0\Library\bin\ffmpeg.exe'

    FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
    ani.save(r'C:\animation.mp4', writer=FFwriter)
else:
    plt.show()






(a,b)=np.polyfit(np.log(t[1:]), np.log(r_max[1:]), deg=1)

if not injection:
    (c,d)=np.polyfit(np.log(t[5000:]), np.log(h_centre[5000:]), deg=1)

plt.loglog(t,h_centre)
if not injection:
    plt.loglog(t[1:],np.exp(d)*t[1:]**c, color='red')

plt.xlabel(r'$t/\tau$')
plt.ylabel('h/H')
plt.title('Hauteur centre')
plt.show()
# plt.plot(t,h_centre)
# plt.xlabel(r'$t/\tau$')
# plt.ylabel('h/H')
# plt.title('Hauteur centre')
# #plt.ylim(0.1,1)
# plt.show()


plt.loglog(t,r_max)

plt.loglog(t[1:],np.exp(b)*t[1:]**a, color='red')
plt.xlabel(r'$t/\tau$')
plt.ylabel(r'$r_{max}/R$')
plt.title('Rayon max')
plt.show()

# plt.plot(np.log(t[1:]),np.log(r_max[1:]))

# plt.title('Rayon max')
# plt.xlabel(r'$t/\tau$')
# plt.ylabel(r'$r_{max}/R$')


plt.show()

plt.plot(t,np.array([volume(h,n) for n in range(n_max)]))
plt.title('Volume total')
if not injection:
    plt.ylim(bottom=0, top=0.00005)
plt.xlabel(r'$t/\tau$')
plt.ylabel(r'$V/V_0$')

plt.show()

theta=np.zeros((n_max,i_max))

for n in range(1,n_max-2):
    for i in range(i_max-1,-1,-1):
        
        iter_theta(theta, n, i)

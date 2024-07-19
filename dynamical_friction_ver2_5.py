# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 17:33:54 2024

@author: Amanda Newton
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

#constants
duration = 1.e8 #years
G = 4 * np.pi**2  #(au^3 / M_sun / yr^2)
M_smbh = 4.e6  #solar masses
m_imbh = 100.  #solar masses
b_max = 1.*206265.  # max deflection parameter, size of cluster
b_min = 0.1 #deflection parameter below, strong deflection
alpha = 1.75  # slope of density profile
rho_0 = 1.35e6  # density in solar masses/pc^3
a0_rho = 0.25 # in pc

#ICs
a0 = 0.01*206265. #au
v0 = np.sqrt((G*M_smbh)/a0) #au/year
L0 = m_imbh*a0*np.sqrt((G*m_imbh)/a0)
dt0 = 2 * np.pi * np.sqrt((a0**3)/(G * M_smbh))
time_elapsed = 0. #years

#initialize lists
dt_list = []
sigma_list = []
ln_Lambd_list = []
rho_list = []
rho_au_list = []
x_list = []
dv_list = []
F_list = []
dL_list = []
L_list = []
a_list = []
v_list = []
t_list = []


#first time for loop
a = a0
a_pc = a/206265.
v = v0
L = L0

#right now, dL is too small, so a_new is too small, so dt gets too small
while time_elapsed < duration and a>1.:
    
    
    #calculate dt, orbital period based on a and M_smbh, and append to dt_list
    #torbital = 2 * np.pi * np.sqrt((a**3)/(G * M_smbh))
    dt = 2 * np.pi * np.sqrt((a**3)/(G * M_smbh))
    #dt = torbital/10.
    dt_list.append(dt)
    # print(f"dt is {dt}")
    # print(f"a is {a}")
    # print(f"v is {v}")
    # print(f"L is {L}")
    
    #functions in dv, dependant on constants and a and v
    sigma = np.sqrt((G*M_smbh)/(a*(1+alpha)))
    ln_Lambd = np.log((a*(sigma**2))/(G*2*m_imbh))
    rho = rho_0 * ((a_pc/a0_rho)**(-alpha))
    rho_au = rho/(206265.**3.)
    x = v/(np.sqrt(2)*sigma)
    
    sigma_list.append(sigma)
    ln_Lambd_list.append(ln_Lambd)
    rho_list.append(rho)
    rho_au_list.append(rho_au)
    x_list.append(x)
    
    #define F, dependant on constants functions and v
    term1 = ((4*np.pi*ln_Lambd*(G**2)*(m_imbh**2)*rho_au)/(v**2))
    term2 = (erf(x) - ((2*x)/np.sqrt(np.pi))*np.exp(-x**2))
    F = term1 * term2
    F_list.append(F)
    
    #forward euler dL, where dL/dt = -Fr, so dependant on F and a
    dL = -F * a * dt
    dL_list.append(dL)
    L_new = L + dL
    L_list.append(L_new)
    # print(f"L_new is {L_new}")
    
    # #solve for a_new from L_new
    # #a_new = (L_new)**2 / (G**2 * m_imbh**3) ### old
    # a_new = (L_new**2.)/(G*M_smbh*(m_imbh**2.))
    # a_list.append(a_new)
    # # print(f"a_new is {a_new}")
    
    # da/a = 2*dL/L
    da = 2.*(dL/L)*a
    a_new = a+da
    a_list.append(a_new)
    
    # #solve for v_new from a_new
    # #v_new = np.sqrt((G*m_imbh)/a_new) ### old
    # v_new = np.sqrt((G*M_smbh)/a_new)
    # v_list.append(v_new)
    # # print(f"v_new is {v_new}")
    
    # v = sqrt(G*M_smbh/a)
    dv = -(da/a)*v/2.
    v_new = v+dv
    v_list.append(v_new)
    
    
    #add dt to t
    time_elapsed += dt
    t_list.append(time_elapsed)
    
    #reset L, a, and v for next loop
    L = L_new
    a = a_new
    #print(f"a is {a}, a_new is {a_new}")
    v = v_new
    
##### adding ICs to lists for plotting
t_list.insert(0,0)
v_list.insert(0,v0)
a_list.insert(0,a0)
dt_list.insert(0,dt0)
L_list.insert(0,L0)
    
plt.scatter(t_list,v_list,s=3)
plt.xlabel("Time [yr]")
plt.ylabel("Velocity [au/year]")
plt.savefig("DF_V_vs_t.png")
plt.show()
plt.scatter(t_list,a_list,s=3)
plt.xlabel("Time [yr]")
plt.ylabel("Semi-Major Axis [au]")
plt.savefig("DF_sma_vs_t.png")
plt.show()
plt.scatter(t_list,dt_list,s=3)
plt.xlabel("Time [yr]")
plt.ylabel("Orbital Period [yr]")
plt.savefig("DF_period_vs_t.png")
plt.show()
plt.scatter(t_list,L_list,s=3)
plt.xlabel("Time [yr]")
plt.ylabel("Angular Momentum [$\mathrm{M}_\odot*au^2/yr$]")
plt.savefig("DF_L_vs_t.png")
plt.show()


# timescale inspiral
T_inspiral = time_elapsed

# timescale relaxation
T_relaxation = (2.e12)/np.log(b_max/b_min)

if T_inspiral == T_relaxation*(1/m_imbh):
    print(f"Yes, T_inspiral is {T_inspiral} and T_relaxation*(1/m_imbh) is {T_relaxation*(1/m_imbh)}")
else:
    print(f"No, T_inspiral is {T_inspiral} and T_relaxation*(1/m_imbh) is {T_relaxation*(1/m_imbh)}")

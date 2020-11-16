# -*- coding: utf-8 -*-
### loop de différent géométrie du glacier

Ltot_list = [5000,10000,20000,40000,80000]
dl = 500
Np_list = [ int(Ltot/dl) for Ltot in Ltot_list ]
results_filename_list = ["ResultsTestdl500Ltot"+str(int(Ltot/1000))+"kmHHT.npz" for Ltot in Ltot_list ]

## reprise du code example_of_simulation_glacial_earthquake.py

import numpy as np
import os
import matplotlib.pyplot as plt
import time
import SeismeGlacier
import AffichageGlacier

vypath = "/home/vyastreb/PEOPLE/Paul_BEGUIN/GIT/StageIPGP/tmp"
PC = "VY" 

if PC != "VY":
    work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modèle SeismeGlacier1'
    results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz\\LoopLtotHHTdl500"
else:
    work_path = vypath 
    results_path = vypath+"/results" 

os.chdir(work_path)

# work space and file name for Fc values

Fc_filename = 'Fc_clear'

# time array
Ttot = 300 
dt = 1
Nt = int(Ttot/dt)
T = np.linspace(0,Ttot,Nt)

# force contact array
Fc = -SeismeGlacier.Fc_lecture(Fc_filename,work_path,dt,Ttot)

# parameters
rho_glace = 920
rho_eau = 1020
alpha_glacier = 1*np.pi/180
g = 9.81
E = 9.3e+9
Cw = 5.623e+6
m = 1/3
H = 800
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_glacier)
ud_eq0 = (rho_glace*H*g*np.sin(alpha_glacier)/Cw)**(1/m)
ud_reg = 0.1*ud_eq0

# For computation process
Nt_r = 100

Ltot_ind = 0
U_bloc_plot = []
Ud_bloc_plot = []
Fsis_bloc_plot = []

for Ltot in Ltot_list :
    
    Np = Np_list[Ltot_ind]
    
    # Appel class Glacier
    glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ud_reg,Np,Ltot,h_im)
    print("Modélisation pour un glacier de longueur Ltot = " + str(int(Ltot/1000)) + 'km' )
    
    # For computation process
    Nt_r = 100 
    
    # mesure  du temps de début de calcul
    t_start = time.time()
    
    # start of the simulation
    alphaHHT = -0.2
    epsHHT = 0.001*ud_eq0
    Ut,Utd,Utdd,Ftf,Ftsismique,Ftsismique_map,Niter_implicite = SeismeGlacier.simu_alpha(glacier,Fc,dt,Nt,alphaHHT,epsHHT,Nt_r)
    
    # end of the computation time
    t_end = time.time()
    print("Temps de calcul: " + "%.2f" % (t_end - t_start))
    
    # saving of the results
    np.savez(results_path + '\\' + results_filename_list[Ltot_ind],Ut=Ut,Utd=Utd,Utdd=Utdd,Ftf=Ftf,Ftsismique=Ftsismique,Ftsismique_map=Ftsismique_map,Niter_implicite=Niter_implicite)
    
    Ltot_ind += 1
    
    U_bloc1_ind = []
    Ud_bloc1_ind = []
    Fsis_bloc1_ind = []
    for k in range(len(Ut)):
        U_bloc1_ind.append(Ut[k][0])
        Ud_bloc1_ind.append(Utd[k][0])
        Fsis_bloc1_ind.append(Ftsismique_map[k][0])
    
    U_bloc_plot.append(U_bloc1_ind)
    Ud_bloc_plot.append(Ud_bloc1_ind)
    Fsis_bloc_plot.append(Fsis_bloc1_ind)

# ======================== #
#  Results for first bloc  #
# ======================== #

plt.figure(1)
for k in range(len(Ltot_list)):
    plt.plot(T,U_bloc_plot[k])

plt.figure(2)
for k in range(len(Ltot_list)):
    plt.plot(T,Ud_bloc_plot[k])

plt.figure(3)
for k in range(len(Ltot_list)):
    plt.plot(T,Fsis_bloc_plot[k])

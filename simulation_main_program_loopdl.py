# -*- coding: utf-8 -*-
### loop de différent géométrie du glacier

Ltot = 5000
dl_list = [1000,500,250,100,50,20]
Np_list = [ int(Ltot/dl) for dl in dl_list ]
results_filename_list = ["ResultsTestdl" + str(dl) + "Ltot5kmHHT.npz" for dl in dl_list ]
results_file = "ResultsTestLoopdlLtot5kmHHT.npz"

## reprise du code example_of_simulation_glacial_earthquake.py

import numpy as np
import os
import matplotlib.pyplot as plt
import time

vypath = "/home/vyastreb/PEOPLE/Paul_BEGUIN/GIT/StageIPGP/tmp"
PC = "VY" 

if PC != "VY":
    work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modèle SeismeGlacier1'
    results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz\\LoopdlHHTLtot5km"
    save_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python'
else:
    work_path = vypath 
    results_path = vypath+"/results" 
    save_path = vypath+"/results" 


os.chdir(work_path)
import SeismeGlacier
import AffichageGlacier

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

for dl in dl_list :
    
    Np = int(Ltot/dl)
    
    # Appel class Glacier
    glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ud_reg,Np,Ltot,h_im)
    print("Modélisation pour une longueur de bloc dl = " + str(dl) + 'm' )
    
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

np.savez(results_path + '\\' + results_file,U=U_bloc_plot,Ud=Ud_bloc_plot,Fsis=Fsis_bloc_plot)


# ======================== #
#  Results for first bloc  #
# ======================== #

file_Fcname = 'Fc_notes.txt'
save_folder = '30 juin Test loop dl'

affichage = AffichageGlacier.Plot_figure_dl_loop(file_Fcname,work_path,results_path,results_file,save_path,save_folder,dl_list,Ttot,ud_reg)
affichage.plot_displacement()
affichage.plot_speed()
affichage.plot_strains()
affichage.plot_error()

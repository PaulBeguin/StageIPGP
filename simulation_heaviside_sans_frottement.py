# -*- coding: utf-8 -*-
### example de simulation pour retournement iceberg et contact du glacier avec le sol


# ========= #
#  Moduls   #
# ========= #

import numpy as np
import os
import matplotlib.pyplot as plt
import time
import os.path

vypath = "/home/vyastreb/PEOPLE/Paul_BEGUIN/GIT/StageIPGP/tmp"
PC = "VY" 

if PC != "VY":
	work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modèle SeismeGlacier sans frottement'
	results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz"
	save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python'
else:
	work_path = vypath
	results_path = os.path.join(vypath,"results")
	save_figure_path = os.path.join(vypath,"results")

os.chdir(work_path)
import SeismeGlacierSansFrott as SeismeGlacier
import AffichageGlacierSansFrott as AffichageGlacier


# ========================= #
#  Loading files T and Fc   #
# ========================= #

# work space and file name for Fc values
Fc_filename = 'Fc_clear'

# time array
Ttot = 300 
dt = 1
Nt = int(Ttot/dt)
T = np.linspace(0,Ttot,Nt)

# force contact array
Fc = np.ones(len(T))
Fc[0] = 1/2 

print("Durée du modèle simulé: " + str(Ttot) + "\n")
print("Nombre de pas de temps: " + str(Nt) + "\n")
print("Pas de temps:" + str(dt) + "\n")


# ================== #
#  Model parameters  #
# ================== #

# Physical
rho_glace = 920
rho_eau = 1020
alpha_ground = 0*np.pi/180 #glacier whithout slope neither at its bedrock nor at its surface
alpha_surface = 0*np.pi/180
g = 9.81
E = 9.3e+9

# Geometric
H = 800
Ltot = 2000
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Np = 20


# Appel class Glacier
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_ground,alpha_surface,g,E,Np,Ltot,h_im)

# For computation process
theme_name = "Heaviside"
results_filename = SeismeGlacier.create_file_result(glacier,results_path,theme_name)
print("Enregistrement des résultats dans fichier " + results_filename + "\n")
Nt_r = 100 


# ======================================= #
#  Simulation et intégration temporelle   #
# ======================================= #

# mesure  du temps de début de calcul
t_start = time.time()

# start of the simulation
alphaHHT = -0.2
epsHHT = 1e-6
Ut,Utd,Utdd,Niter_implicite = SeismeGlacier.simu_alpha(glacier,Fc,dt,Nt,alphaHHT,epsHHT,Nt_r)

# end of the computation time
t_end = time.time()
print("Temps de calcul: " + "%.2f" % (t_end - t_start))

# saving of the results
np.savez(os.path.join(results_path,results_filename),Ut=Ut,Utd=Utd,Utdd=Utdd,Niter_implicite=Niter_implicite)
	


# ========================= #
#  Affichage de résultats   #
# ========================= #

# récupération de Nt_retour et T_retour
save_folder_name = SeismeGlacier.create_folder_figure(glacier,save_figure_path,theme_name)
print("Enregistrement des figures dans fichier " + save_folder_name + "\n")

lu_plot=[0,500,1000,1500,Ltot]

affichage = AffichageGlacier.Plot_figures(results_filename,work_path,results_path,save_figure_path,save_folder_name,lu_plot,Ttot,Ltot,Fc)
affichage.plot_displacement()
affichage.plot_speed()
# affichage.plot_strains()
affichage.plot_Niter()
affichage.get_map()
affichage.map_displacement()
# affichage.map_displacement2()
# affichage.map_sismique()

### example de simulation pour retournement iceberg et contact du glacier avec le sol


# ========= #
#  Moduls   #
# ========= #

import numpy as np
import os
import matplotlib.pyplot as plt
import time

work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modèle SeismeGlacier1'
os.chdir(work_path)
import SeismeGlacier
import AffichageGlacier


# ========================= #
#  Loading files T and Fc   #
# ========================= #

# work space and file name for Fc values
results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz"
Fc_filename = 'Fc_clear'

# time array
Ttot = 300 
dt = 1
Nt = int(Ttot/dt)
T = np.linspace(0,Ttot,Nt)

# force contact array
Fc = -SeismeGlacier.Fc_lecture(Fc_filename,work_path,dt,Ttot)

print("Durée du modèle simulé: " + str(Ttot) + "\n")
print("Nombre de pas de temps: " + str(Nt) + "\n")
print("Pas de temps:" + str(dt) + "\n")


# ================== #
#  Model parameters  #
# ================== #

# Physical
rho_glace = 920
rho_eau = 1020
alpha_ground = 1*np.pi/180
alpha_surface = 1*np.pi/180
g = 9.81
E = 9.3e+9
Cw = 5.623e+6
m = 1/3
Cc = 0.4
type_law = 0 #0=Weertman, 1=Tsai 

# Geometric
H = 800
Ltot = 2000
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Np = 20

# Steady sliding
ud_eq0 = (rho_glace*H*g*np.sin(alpha_ground)/Cw)**(1/m)
ud_reg = 0.1*ud_eq0

# Appel class Glacier
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_ground,alpha_surface,g,E,Cw,m,Cc,type_law,ud_reg,Np,Ltot,h_im)

# For computation process
theme_name = "Test"
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
epsHHT = 0.001*ud_eq0
Ut,Utd,Utdd,Ftf,Ftsismique,Ftsismique_map,Niter_implicite = SeismeGlacier.simu_alpha(glacier,Fc,dt,Nt,alphaHHT,epsHHT,Nt_r)

# end of the computation time
t_end = time.time()
print("Temps de calcul: " + "%.2f" % (t_end - t_start))

# saving of the results
np.savez(results_path + '\\' + results_filename,Ut=Ut,Utd=Utd,Utdd=Utdd,Ftf=Ftf,Ftsismique=Ftsismique,Ftsismique_map=Ftsismique_map,Niter_implicite=Niter_implicite)


# ========================= #
#  Affichage de résultats   #
# ========================= #

# récupération de Nt_retour et T_retour
Fc_note_filename = 'Fc_notes.txt'
save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python'
save_folder_name = SeismeGlacier.create_folder_figure(glacier,save_figure_path,theme_name)
print("Enregistrement des figures dans fichier " + save_folder_name + "\n")

lu_plot=[0,500,1000,1500,Ltot]

affichage = AffichageGlacier.Plot_figures(results_filename,Fc_note_filename,work_path,results_path,save_figure_path,save_folder_name,lu_plot,Ttot,Ltot,Fc,ud_reg)
affichage.plot_displacement()
affichage.plot_speed()
affichage.plot_strains()
affichage.plot_Niter()
affichage.map_displacement()
# affichage.map_displacement2()
affichage.map_sismique()

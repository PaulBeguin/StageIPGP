### example de simulation pour retournement iceberg et contact du glacier avec le sol

<<<<<<< HEAD

# ========= #
#  Moduls   #
# ========= #
=======
# ========== #
#  Moduls   #
# ========== #
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe

import numpy as np
import os
import matplotlib.pyplot as plt
import time

work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modèle SeismeGlacier1'
os.chdir(work_path)
import SeismeGlacier
import AffichageGlacier

<<<<<<< HEAD

# ========================= #
#  Loading files T and Fc   #
# ========================= #

# work space and file name for Fc values
results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz"
results_filename = "ResultsTestdl100Ltot2kmHHT.npz"
=======
# =============================== #
#  Loading files T and Fc   #
# =============================== #

# work space and file name for Fc values
results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results Copy File Npz"
results_filename = "ResultsTestdl500Ltot4kmHHT.npz"
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
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

<<<<<<< HEAD

=======
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
# ================== #
#  Model parameters  #
# ================== #

# Physical
rho_glace = 920
rho_eau = 1020
alpha_glacier = 1*np.pi/180
g = 9.81
E = 9.3e+9
Cw = 5.623e+6
m = 1/3

# Geometric
H = 800
<<<<<<< HEAD
Ltot = 2000
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_glacier)
Np = 20

# Steady sliding
ud_eq0 = (rho_glace*H*g*np.sin(alpha_glacier)/Cw)**(1/m)
ud_reg = 0.1*ud_eq0

# Appel class Glacier
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ud_reg,Np,Ltot,h_im)
=======
Ltot = 40000
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_glacier)
Np = 80
Np_r = 1 # à supprimer

# Steady sliding
ud_eq0 = ( rho_glace*H*g*np.sin(alpha_glacier)/Cw)**(1/m)
ud_reg = 0.1*ud_eq0

# Appel class Glacier
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ud_reg,Np,Ltot,h_im,Np_r)
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe

# For computation process
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

<<<<<<< HEAD

=======
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
# ========================= #
#  Affichage de résultats   #
# ========================= #

# récupération de Nt_retour et T_retour
Fc_note_filename = 'Fc_notes.txt'
save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python'
<<<<<<< HEAD
save_folder_name = '30 juin Ltot2km'
lu_plot=[0,500,1000,1500,Ltot]

affichage = AffichageGlacier.Plot_figures(results_filename,Fc_note_filename,work_path,results_path,save_figure_path,save_folder_name,lu_plot,Ttot,Ltot,Fc,ud_reg)
affichage.plot_displacement()
affichage.plot_speed()
affichage.plot_strains()
affichage.plot_Niter()
=======
save_folder_name = '26 juin HHT test'
lu_plot=[0,500,1000,1500,Ltot]

affichage = AffichageGlacier.Plot_figures(results_filename,Fc_note_filename,work_path,results_path,save_figure_path,save_folder_name,lu_plot,Ttot,Ltot,Fc,ud_reg)
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe

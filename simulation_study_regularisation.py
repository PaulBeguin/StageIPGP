### example de simulation pour retournement iceberg et contact du glacier avec le sol


# ========= #
#  Moduls   #
# ========= #

import numpy as np
import os
import matplotlib.pyplot as plt
import time

vypath = "/home/vyastreb/PEOPLE/Paul_BEGUIN/GIT/StageIPGP/tmp"
PC = "PB" 

if PC != "VY":
    work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modele_SeismeGlacier_regularisation'
    results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Results_Copy_File_Npz"
    save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python'
else:
    work_path = vypath 
    results_path = vypath+"/results" 
    save_path = vypath+"/results" 


os.chdir(work_path)
import SeismeGlacierReg as SeismeGlacier
import AffichageGlacierReg as AffichageGlacier


# ========================= #
#  Loading files T and Fc   #
# ========================= #

# work space and file name for Fc values
Fc_filename = 'Fc_clear'

# time array
Ttot = 300 
# T_elas = 2*np.pi*np.sqrt(glacier.rho_glace/glacier.E)*np.min(glacier.dl_l)
# dt = T_elas/10
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
dl = Ltot/Np

# For computation process
theme_name = "LoopUdReg"
results_filename = SeismeGlacier.create_file_result(dl,Ltot,type_law,results_path,theme_name)
print("Enregistrement des résultats dans fichier " + results_filename + "\n")
Nt_r = 100

# Steady sliding
ud_eq0 = (rho_glace*H*g*np.sin(alpha_ground)/Cw)**(1/m)
ud_reg_list = [0.2,0.1,0.05,0.02,0.01,0.005] #list of coeffiecient for set the ud_reg parameter

Ut_result0 = []
Utd_result0 = []
Ft_result0 = []
Niter_result = []

for k in range(len(ud_reg_list)): #loop of Ud_reg 

    ud_reg = ud_reg_list[k]*ud_eq0
    print('## Itération numéro : ' + str(k) + "\n## Vitesse de régularisation : " + "%.3f" % (ud_reg*1e+6) + "mu m/s")

    # Call of Glacier classe in SeismeGlacier
    glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H,alpha_ground,alpha_surface,g,E,Cw,m,Cc,type_law,ud_reg,Np,Ltot,h_im)
    
    
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
    
    Ut_result0k, Utd_result0k, Ft_result0k = [], [], []
    for i in range(len(T)):
        Ut_result0k.append(Ut[i][0])
        Utd_result0k.append(Utd[i][0])
        Ft_result0k.append(Ftsismique_map[i][0])
    
    Ut_result0.append(Ut_result0k)
    Utd_result0.append(Utd_result0k)
    Ft_result0.append(Ft_result0k)
    Niter_result.append(Niter_implicite)
    
# saving of the results
np.savez(results_path + '\\' + results_filename,Ut0=Ut_result0,Utd0=Utd_result0,Ft0=Ft_result0,Niter=Niter_result)


# ========================= #
#  Affichage de résultats   #
# ========================= #

# récupération de Nt_retour et T_retour
Fc_note_filename = 'Fc_notes.txt'

save_folder_name = SeismeGlacier.create_folder_figure(save_figure_path,theme_name)
print("Enregistrement des figures dans fichier " + save_folder_name + "\n")

affichage = AffichageGlacier.Plot_figures(results_filename,Fc_note_filename,work_path,results_path,save_figure_path,save_folder_name,ud_reg_list,Ttot,Ltot,Fc,ud_eq0)
affichage.plot_displacement()
# affichage.plot_speed()
affichage.plot_strains()
affichage.plot_Niter()
# affichage.map_displacement()
# affichage.map_displacement2()
# affichage.map_sismique()

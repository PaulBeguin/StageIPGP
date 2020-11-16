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
    work_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modele_SeismeCoulomb3'
    results_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\2020_11_5_Results_10blocs_Coulomb"
    save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_11_5_CoulWeert_Glacier'
else:
    work_path = vypath 
    results_path = vypath+"/results" 
    save_path = vypath+"/results" 

os.chdir(work_path)
import SeismeGlacier3 as SeismeGlacier


# ================== #
#  Model parameters  #
# ================== #

# Physical
rho_glace = 920
rho_eau = 1025
alpha_ground = 1*np.pi/180
alpha_surface = alpha_ground
g = 9.81
E = 9.3e+9
ud_eq0 = 10/(3600*24)
m = 1/3
Cc = 0.5
type_law = 1 #0=Weertman, 1=Coulomb

# Geometric
H = 800
Ltot = 1400
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Xt = 0
Xc = -1400
Np = 10
H_filename = 'H_glacier'

# Steady sliding
Cw = (rho_glace*H*g*np.sin(alpha_ground))/((ud_eq0)**(m))
ud_reg = 0.1*ud_eq0

# Appel class Glacier
Nh = 2
# writing of the corner coordinates of the glacier and its tongue
SeismeGlacier.Write_coord_file_linear(alpha_ground,alpha_surface,Xt,Ltot,H,h_im,H_filename)
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H_filename,alpha_ground,alpha_surface,g,E,Cw,m,Cc,type_law,ud_reg,Np,Ltot,H,h_im,Xc)


# initial state
Ut0_init = np.zeros(Np+1)
Utd0_init = ud_eq0*np.ones(Np+1)
Utdd0_init = np.zeros(Np+1)

# Filename for results
results_filename_statique = "Results10blocsCoul_statique"
results_filename = "Results10blocsCoul"
print("Enregistrement des résultats dans fichier " + results_filename + "\n")
Nt_r = 1

# ========================= #
#  Loading files T and Fc   #
# ========================= #

# work space and file name for Fc values
Fc_filename = 'Fc_clear'

# time array
M = glacier.M[:-1,:-1]*10000
K = glacier.K[:-1,:-1]

Ttot = 150
T_start = 20
T_retour = 84.84
dt = 0.001
Nt = int(Ttot/dt)+1
T = np.linspace(0,Ttot,Nt)
Ut_last = ud_eq0*T

# force contact array
Fc0 = np.zeros(Nt)
Fc = SeismeGlacier.Fc_lecture(Fc_filename,work_path,T,T_start)

print("Durée du modèle simulé: " + str(Ttot) + "\n")
print("Nombre de pas de temps: " + str(Nt) + "\n")
print("Le pas de temps est de " + "%.6f" % dt + ", régulier \n")

# ======================================= #
#  Simulation et intégration temporelle   #
# ======================================= #

# start of the simulations
alpha_HHT = -1/3
eps_HHT = g/1000

# static
print("Début calcul intégration en temps")
t_start = time.time()

# définition des vecteurs de calcul et de l'état initial du système
Un = Ut0_init # displacement
Udn = Utd0_init # vitesse 
Uddn = Utdd0_init # acceleration

# arrays to return and at plotting 
Ut = [[Un[k] for k in range(Np+1)]] # displacement
Utd = [[Udn[k] for k in range(Np+1)]] # vitesse affichée
Utdd = [[Uddn[k] for k in range(Np+1)]] # acceleration

Fttot_map = [[0 for k in range(Np+1)]] # force de frottement comparée au frottement statique
Ftsismique_map = [[0 for k in range(Np+1)]]
Adh_first = [False for k in range(Np+1)]
sign_last = [ 1 for k in range(Np+1)]

Niter_implicite = []

# matrix and coefficient
mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

#calcul of the first strain vectors
Poidsx = -glacier.P_xpoids[:-1]

Ff_coul = np.zeros(Np+1)
Ff_tot = np.zeros(Np+1)

Adh = [False for k in range(Np+1)]
sign_last = np.ones(Np+1)

# force de perturbation
Fp = np.zeros(Np+1)
Fp[0] = Fc[0]
# force élastique avec commposante au dernier bloc - bloc rigide qui se déplace à la vitesse ud_eq0
Fe = - np.dot(K,Un)
Fe_last = np.zeros(Np+1)
Fe_last[-1] = -K[-1][-1]*(Un[-1] - Ut_last[0])
# force de frottement nulle 
Ff = np.zeros(Np+1)

P_n = Poidsx + Fe + Fe_last + Fp + Ff

for n in range(Nt-1):
    
    un_i = Un + dt*Udn + 0.5*(dt**2)*Uddn
    udn_i = Udn + dt*Uddn
    uddn_i = Uddn
    
    Ms_i = M + mk*K
    
    Fe_i = - np.dot(K,un_i)
    Fe_last_i = Fe_last.copy()
    Fe_last_i[-1] = -K[-1][-1]*(un_i[-1] - Ut_last[n+1])
    Fp_i = Fp.copy()
    Fp_i[0] = Fc[n+1]
    
    F_trial_i = Poidsx + Fe_i + Fe_last_i + Fp_i
    
    for k in range(Np+1):
        if Adh[k] == False :
            if udn_i[k]==0 or Udn[k]*udn_i[k]<0: # passage de glissement à adhérence
                # force de frottement à définir en supposant l'état d'adhérence
                Adh[k] = True
        elif Adh[k] == True :
            if np.abs(F_trial_i[k]) > glacier.F_coulomb_lim[k] : # passage d'adhérence à glissement
                Adh[k] = False
                if un_i[k] != 0:
                    sign_last[k] = np.sign(un_i[k])
                else :
                    sign_last[k] = np.sign(F_trial_i[k])
    
    Ff_alpha = np.zeros(Np+1)
    P_alpha_t = np.zeros(Np+1)
    P_alpha_t = F_trial_i
    # for k in range(Np+1):
    #     if Adh[k] == True : # sticking state at k bloc
    #         Ff_alpha[k] = -F_trial_i[k]*sign_last[k]
    #         P_alpha_t[k] = 0
    #     if Adh[k] == False : # sliping state at k bloc
    #         Ff_alpha[k] = -glacier.F_coulomb_lim[k]*sign_last[k]
    #         P_alpha_t[k] = Ff_alpha[k] + F_trial_i[k]
    
    r_i = p_alpha*(P_n) + p_alpha1*(P_alpha_t) - np.dot(M,uddn_i)
    delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
    
    un_i1 = un_i + uiai*delta_uddn_i
    udn_i1 = udn_i + udiai*delta_uddn_i
    uddn_i1 = uddn_i + delta_uddn_i
    for k in range(Np+1):
        if Adh[k] == True :
            un_i1[k] = un_i[k]
            udn_i1[k] = 0
            uddn_i1[k] = 0
    
    i = 0
    
    while i<50 and np.max(np.abs(uddn_i1 - uddn_i))>eps_HHT :
        
        Fe_i = - np.dot(K,un_i1)
        Fe_last_i = Fe_last.copy()
        Fe_last_i[-1] = -K[-1][-1]*(un_i1[-1] - Ut_last[n+1])
        Fp_i = Fp.copy()
        Fp_i[0] = Fc[n+1]
        
        F_trial_i = Poidsx + Fe_i + Fe_last_i + Fp_i
        
        for k in range(Np+1):
            if Adh[k] == False :
                if udn_i1[k]==0 or Udn[k]*udn_i1[k]<0: # passage de glissement à adhérence
                    # force de frottement à définir en supposant l'état d'adhérence
                    Adh[k] = True
            elif Adh[k] == True :
                if np.abs(F_trial_i[k]) > glacier.F_coulomb_lim[k] : # passage d'adhérence à glissement
                    Adh[k] = False
                    if un_i1[k] != 0:
                        sign_last[k] = np.sign(un_i1[k])
                    else :
                        sign_last[k] = np.sign(F_trial_i[k])
        
        Ff_alpha = np.zeros(Np+1)
        P_alpha_t = np.zeros(Np+1)
        P_alpha_t = F_trial_i
        # for k in range(Np+1):
        #     if Adh[k] == True : # sticking state at k bloc
        #         Ff_alpha[k] = -F_trial_i[k]*sign_last[k]
        #         P_alpha_t[k] = 0
        #     if Adh[k] == False : # sliping state at k bloc
        #         Ff_alpha[k] = -glacier.F_coulomb_lim[k]*sign_last[k]
        #         P_alpha_t[k] = Ff_alpha[k] + F_trial_i[k]
        
        un_i = un_i1
        udn_i = udn_i1
        uddn_i = uddn_i1
        
        Ms_i = M + mk*K
        
        r_i = p_alpha*(P_n) + p_alpha1*(P_alpha_t) - np.dot(M,uddn_i)
        delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
        
        un_i1 = un_i + uiai*delta_uddn_i
        udn_i1 = udn_i + udiai*delta_uddn_i
        uddn_i1 = uddn_i + delta_uddn_i
        for k in range(glacier.lim_ind_coul,glacier.lim_ind_weert+1):
            if Adh[k] == True :
                un_i1[k] = un_i[k]
                udn_i1[k] = 0
                uddn_i1[k] = 0
        
        i +=1
    
    en = np.max(np.abs(udn_i1 - udn_i)) #residual error
    
    # state vector at n+1
    Un = un_i1
    Udn = udn_i1
    Uddn = uddn_i1
    
    P_n = P_alpha_t
    Ff_n = Ff_alpha
    
    # enregistrement des données pour affichage
    if n%1000 == 0:
        print(n) 
    
    if i==50 :
        print('Alerte non convergence au pas ' + str(n))
    
    Ut.append( Un.tolist() )
    Utd.append( Udn.tolist() )
    Utdd.append( Uddn.tolist() )
    
    Ftsismique_map.append( Ff_n.tolist() )
    Fttot_map.append( P_n.tolist() )
    Niter_implicite.append(i)
    

t_end = time.time() 
print("Temps de calcul statique: " + "%.2f" % (t_end - t_start) + '\n')

# T_r = np.linspace(0,Ttot,len(Ut))
# np.savez(results_path + '\\' + results_filename + str(ind_memo) + 'T_r.npz',T_r=T_r)
np.savez(results_path + '\\' + results_filename + 'Ut.npz',Ut=np.array(Ut))
np.savez(results_path + '\\' + results_filename_statique + 'Utd.npz',Utd=np.array(Utd))
# np.savez(results_path + '\\' + results_filename_statique + str(ind_memo) + 'Utdd.npz',Utdd=np.array(Utdd))
np.savez(results_path + '\\' + results_filename_statique + 'Ftsismique_map.npz',Ftsismique_map=np.array(Ftsismique_map))
# np.savez(results_path + '\\' + results_filename_statique + str(ind_memo) + 'Ft_map.npz',Ft_map=np.array(Fttot_map))
np.savez(results_path + '\\' + results_filename_statique + 'Niter_implicite.npz',Niter_implicite=np.array(Niter_implicite))


# del T_r,Ut,Utd,Utdd,Ftsismique_map,Fttot_map,Niter_implicite








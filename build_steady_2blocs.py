# test pour un bloc 

import numpy as np
import matplotlib.pyplot as plt
import os

work_path = "C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modele_SeismeCoulomb3"
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
Ud0_last= -10/(3600*24)
m = 1/3
Cc = np.tan(alpha_ground)
type_law = 1 #0=Weertman, 1=Coulomb

# Geometric
H = 800
Ltot = 200/np.cos(alpha_ground)
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Xt = -1300
Xc = -1400
Np = 1
H_filename = 'H_glacier'

# Steady sliding
Cw = (rho_glace*H*g*np.sin(alpha_ground))/((np.abs(Ud0_last))**(m))
ud_reg = 0.1*np.abs(Ud0_last)

# Appel class Glacier
SeismeGlacier.Write_coord_file_linear(alpha_ground,alpha_surface,Xt,Ltot,H,h_im,H_filename)
glacier = SeismeGlacier.Glacier(rho_glace,rho_eau,H_filename,alpha_ground,alpha_surface,g,E,Cw,m,Cc,type_law,ud_reg,Np,Ltot,H,h_im,Xc)

X_down = np.array([Xt,Xc,Xt-Ltot*np.cos(alpha_ground)])
dl_l = np.array([X_down[0]-X_down[1],X_down[1]-X_down[2]])

M_bloc = rho_glace*H*Ltot/2
Poids_x = -np.sin(alpha_ground)*M_bloc*g* np.array([1,1])
Poids_y = -np.cos(alpha_ground)*M_bloc*g* np.array([1,1])
# P_gli = glacier.F_coulomb_lim[0]
P_gli = Cc*np.abs(Poids_y[0])

eps_HHT = P_gli/(100*M_bloc)

K = E*H/(Ltot/2)*np.array([[1,-1],[-1,1]])
M = M_bloc*np.array([[1,0],[0,1]])
T_end = 30
T_Fp = 15
dt = 0.01

Nt = int(T_end/dt) +1
T = np.linspace(0,T_end,Nt)
Fp = np.zeros(Nt)

# lancement boucle temporelle
alpha_HHT = -0.005
U0,Ud0,Udd0 = np.array([0,0]), np.array([0,0]), np.array([0,0])
U,Ud,Udd = U0,Ud0,Udd0
Ff = np.array([0,0])
U_plot,Ud_plot,Udd_plot = [U0],[Ud0],[Udd0]
Ff_plot = []
F_plot = []
F_trial_plot = []
Niter_implicite,ErrResi = [],[]

Nt = int(T_end/dt)+1
T = np.linspace(0,T_end,Nt)
U_last = Ud0_last*T
Fe_last = np.zeros(2)
Fp = np.zeros(Nt)
p = 0
while p<Nt :
    if T[p] >= T_Fp :
        Fp[p] = P_gli*np.sin(np.pi*(T[p]-T_Fp))
    p+=1

# matrix and coefficient
mc = (1+alpha_HHT)*(0.5-alpha_HHT)*dt
mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

Adh = [False,False]
Adh_first = [True,True]
Adh_list = []
ind_zoom = 0
sign_last = np.array([1,1])

Fe_last[1] = -K[1][1]*(U[1] - U_last[0])
P_n = -np.dot(K,U) + Poids_x + Fp[0] + Ff + Fe_last

for n in range(Nt-1):
    
    # print("PAS :" +str(n))
    
    dt = T[n+1] - T[n] #variable step time
    Un = U
    Udn = Ud
    Uddn = Udd
    
    U_i = U + dt*Ud + 0.5*(dt**2)*Udd
    Ud_i = Ud + dt*Udd
    Udd_i = Udd
    
    C_i = np.zeros((2,2))
    if np.abs(Ud_i[1]) > ud_reg :
        C_i[1][1] = glacier.Cw*glacier.m*dl_l[1]*(np.abs(Ud_i[1])**(glacier.m-1))
    else :
        C_i[1][1] = glacier.Cl*dl_l[1]
    
    Ms_i = M + mc*C_i + mk*K
    
    ## calcul de glacier.Ff_coul2(un_i,udn_i,Udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
    Fe_last[1] = -K[1][1]*(U[1] - U_last[n+1])
    F_trial_i = -np.dot(K,np.transpose(U_i)) + Poids_x + Fp[n+1] + Fe_last
    
    ## définition des états de glissement et adhérence au premier bloc de Coulomb
    for k in range(Np+1):
        if (X_down[k]<=0) and (X_down[k]>Xc):
            if Adh[k] == False :
                if Ud_i[k]==0 or Udn[k]*Ud_i[k]<0: # passage de glissement à adhérence
                    # force de frottement à définir en supposant l'état d'adhérence
                    Adh[k] = True
            elif Adh[k] == True :
                if np.abs(F_trial_i[k]) > P_gli : # passage d'adhérence à glissement
                    Adh[k] = False
                    if Ud_i[0] != 0:
                        sign_last[k] = np.sign(Ud_i[k])
                    else :
                        sign_last[k] = np.sign(F_trial_i[k])
    # print(Adh)
    
    ## frottement de Weertman pour le deuxième bloc
    Ff_weert = np.array([0,0])
    for k in range(Np+1):
        if X_down[k] <= Xc :
            if np.abs(Ud_i[k])>ud_reg :
                Ff_weert[k]=-np.sign(Ud_i[k])*glacier.Cw*dl_l[k]*(np.abs(Ud_i[k])**(glacier.m))
            else :
                Ff_weert[k]=-glacier.Cl*dl_l[k]*Ud_i[k]
    
    ## calcul des forces Coulomb et calcul des résultantes
    Ff_alpha = Ff_weert # initialisation avec le frottement aux blocs de Weertman seulement
    P_nalpha_t = np.zeros(glacier.Np-1)
    for k in range(glacier.lim_ind_coul,glacier.lim_ind_weert):
        if Adh[k] == True : #adhérence , la résultante au noeud devient nulle et la force de frottement est compense les autres efforts au noeud
            Ff_alpha[k] = -F_trial_i[k]*sign_last[k]
            P_nalpha_t[k] = 0
        if Adh[k] == False : #glissement, la résultante au noeud égale à la somme des efforts extérieur et de frottement, la force de frottement reste égale à celle à l'état précédent
            Ff_alpha[k] = -P_gli*sign_last[k]
            P_nalpha_t[k] = Ff_alpha[k] + F_trial_i[k]
    
    for k in range(glacier.lim_ind_weert,glacier.Np-1) :
        P_nalpha_t[k] = F_trial_i[k] + Ff_alpha[k]
    
    # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
    r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha_t) - np.dot(M_bloc,Udd_i)
    delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
    
    U_i1 = U_i + uiai*delta_uddn_i
    Ud_i1 = Ud_i + udiai*delta_uddn_i
    Udd_i1 = Udd_i + delta_uddn_i
    
    for k in range(Np+1):
        if Adh[k] == True :
            U_i1[k] = U_i[k]
            Ud_i1[k] = 0
            Udd_i1[k] = 0
    
    i = 0
    
    while i<10 and np.max(np.abs(Udd_i1-Udd_i))>eps_HHT :
        
        # print('Itération :' + str(i+1))
        
        C_i = np.zeros((2,2))
        if np.abs(Ud_i1[1]) > ud_reg :
            C_i[1][1] = glacier.Cw*glacier.m*dl_l[1]*(np.abs(Ud_i1[1])**(glacier.m-1))
        else :
            C_i[1][1] = glacier.Cl*dl_l[1]
        
        Ms_i = M + mc*C_i + mk*K
        
        ## calcul de glacier.Ff_coul2(un_i,udn_i,Udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
        Fe_last[1] = -K[1][1]*(U[1] - U_last[n+1])
        F_trial_i = -np.dot(K,np.transpose(U_i1)) + Poids_x + Fp[n+1] + Fe_last
        
        ## définition des états de glissement et adhérence au premier bloc de Coulomb
        for k in range(Np+1):
            if (X_down[k]<=0) and (X_down[k]>Xc):
                if Adh[k] == False :
                    if Ud_i1[k]==0 or Udn[k]*Ud_i1[k]<0: # passage de glissement à adhérence
                        # force de frottement à définir en supposant l'état d'adhérence
                        Adh[k] = True
                elif Adh[k] == True :
                    if np.abs(F_trial_i[k]) > P_gli : # passage d'adhérence à glissement
                        Adh[k] = False
                        if Ud_i1[0] != 0:
                            sign_last[k] = np.sign(Ud_i1[k])
                        else :
                            sign_last[k] = np.sign(F_trial_i[k])
        
        # print(Adh)
        ## frottement de Weertman pour le deuxième bloc
        Ff_weert = np.array([0,0])
        for k in range(Np+1):
            if X_down[k] <= Xc :
                if np.abs(Ud_i1[k])>ud_reg :
                    Ff_weert[k]=-np.sign(Ud_i1[k])*glacier.Cw*dl_l[k]*(np.abs(Ud_i1[k])**(glacier.m))
                else :
                    Ff_weert[k]=-glacier.Cl*dl_l[k]*Ud_i1[k]
        
        U_i = U_i1
        Ud_i = Ud_i1
        Udd_i = Udd_i1
        
        ## calcul des forces Coulomb et calcul des résultantes
        Ff_alpha = Ff_weert # initialisation avec le frottement aux blocs de Weertman seulement
        P_nalpha_t = np.zeros(glacier.Np-1)
        for k in range(glacier.lim_ind_coul,glacier.lim_ind_weert):
            if Adh[k] == True : #adhérence , la résultante au noeud devient nulle et la force de frottement est compense les autres efforts au noeud
                Ff_alpha[k] = -F_trial_i[k]*sign_last[k]
                P_nalpha_t[k] = 0
            if Adh[k] == False : #glissement, la résultante au noeud égale à la somme des efforts extérieur et de frottement, la force de frottement reste égale à celle à l'état précédent
                Ff_alpha[k] = -P_gli*sign_last[k]
                P_nalpha_t[k] = Ff_alpha[k] + F_trial_i[k]
        
        for k in range(glacier.lim_ind_weert,glacier.Np-1) :
            P_nalpha_t[k] = F_trial_i[k] + Ff_alpha[k]
        
        # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
        r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha_t) - np.dot(M_bloc,Udd_i)
        delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
    
        U_i1 = U_i + uiai*delta_uddn_i
        Ud_i1 = Ud_i + udiai*delta_uddn_i
        Udd_i1 = Udd_i + delta_uddn_i
        
        for k in range(glacier.lim_ind_coul,glacier.lim_ind_weert):
            if Adh[k] == True :
                U_i1[k] = U_i[k]
                Ud_i1[k] = 0
                Udd_i1[k] = 0
        
        i +=1
    
    if i==10 :
        print('Alerte non convergence au pas ' + str(n))
    
    Ff = Ff_alpha
    P_n = P_nalpha_t
    
    en = np.abs(Ud_i1 - Ud_i) #residual error
    
    # state vector at n+1
    U = U_i1
    Ud = Ud_i1
    Udd = Udd_i1
    
    U_plot.append( U )
    Ud_plot.append( Ud )
    Udd_plot.append( Udd )
    
    Ff_plot.append( Ff )
    F_plot.append(P_nalpha_t)
    F_trial_plot.append(F_trial_i)
    Adh_list.append( Adh )
    
    Niter_implicite.append(i)
    ErrResi.append(en)

save_figure_path = 'C:\\Users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_10_19_CoulWeert_2blocs_test_CcTan'
os.chdir(save_figure_path)

import matplotlib
matplotlib.rcParams['figure.figsize'] = 35,18
matplotlib.rcParams.update({'font.size': 20})
title_font = {'size':'30'}
label_font = {'size':'30'}

Cc_label = 'Tan'
plt.figure(1)
plt.plot(T,np.array(U_plot)*1e+3,linewidth=2)
plt.plot(T,U_last*1e+3,'b--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Déplacement ($mm$)',label_font)
plt.title("Déplacement du modèle 2 blocs avec $\mu = tan( alpha ) $",label_font)
plt.grid(True)
# plt.axis([0,20,-5,5])
plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0 .T$'],loc='best')
plt.savefig('DeplacementTestCc' + Cc_label + '.png')
plt.savefig('DeplacementTestCc' + Cc_label + '.pdf')

plt.figure(2)
plt.plot(T,np.array(Ud_plot)*1e+3,linewidth=2)
plt.plot(np.array([T[0],T[-1]]),np.array([Ud0_last*1e+3,Ud0_last*1e+3]),'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Vitesse ($mm/s$)',label_font)
plt.title("Vitesse du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0$'],loc='best')
plt.grid(True)
# plt.axis([0,20,-5,5])
plt.savefig('VitesseTestCc' + Cc_label + '.png')
plt.savefig('VitesseTestCc' + Cc_label + '.pdf')

plt.figure(3)
plt.plot(T[1:],np.array(Ff_plot)*1e-6,linewidth=2)
plt.plot(T,Fp*1e-6,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.plot(np.array([T[0],T[-1]]),np.array([-P_gli*1e-6,-P_gli*1e-6]),'b--',linewidth=2)
plt.plot(np.array([T[0],T[-1]]),np.array([P_gli*1e-6,P_gli*1e-6]),'b--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Force de frottement ($MN/m$)',label_font)
plt.title("Frottement du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.legend(['Bloc Coulomb','Bloc Weertman','Force de perturbation','Limites $-P_{gli}$ et $P_{gli}$ de glissement'],loc='best')
plt.grid(True)
plt.savefig('FrottementTestCc' + Cc_label + '.png')
plt.savefig('FrottementTestCc' + Cc_label + '.pdf')

plt.figure(4)
plt.plot(T[1:],np.array(F_plot)*1e-6,linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Resultante ($MN/m$)',label_font)
plt.title("Résultante du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.grid(True)
# plt.axis([0,20,-10,10])
plt.legend(['Bloc Coulomb','Bloc Weertman'],loc='best')
plt.savefig('ResultanteTestCc' + Cc_label + '.png')
plt.savefig('ResultanteTestCc' + Cc_label + '.pdf')

# test pour un bloc 

import numpy as np
import matplotlib.pyplot as plt
import os

work_path = "C:\\users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Programmes Python\\Modele_SeismeCoulomb3"
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
ud0_last= -10/(3600*24)
m = 1/3
Cc = np.tan(alpha_ground)
type_law = 1 #0=Weertman, 1=Coulomb

# Geometric
H = 800
Ltot = 1000/np.cos(alpha_ground)
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Xt = -400
Np = 2

# goemetric definition for the glaciar
X_down = np.linspace(Xt,Xt-Ltot,Np+1)
dl_l = np.array([X_down[k]-X_down[k+1] for k in range(len(X_down)-1)])
dl_liai = np.array([(dl_l[k]+dl_l[k+1])/2 for k in range(len(dl_l)-1)])
dl_min = np.min(dl_l)

# mass and weight for the blocs
M_bloc = rho_glace*H*dl_l
Poids_x = -np.sin(alpha_ground)*g*M_bloc
Poids_y = -np.cos(alpha_ground)*g*M_bloc
P_gli = Cc*np.abs(Poids_y)

# convergence criteria
eps_HHT = np.min(P_gli/(100*M_bloc))

# stiffness matrix
K = np.zeros((Np,Np))
ke = E*H*np.array([[1,-1],[-1,1]])
for i in range(Np-1):
    K[i:i+2,i:i+2] = K[i:i+2,i:i+2] + ke/dl_liai[i]

# mass matrix
M = np.diag(M_bloc)

# time simulated - discretisation
T_end = 5
T_Fp = 2
dt = 0.1*dl_min*np.sqrt(rho_glace/E)
Nt = int(T_end/dt) + 1
T = np.linspace(0,T_end,Nt)

# definition des variables du problème
u0,ud0,udd0 = np.zeros(Np), ud0_last*np.ones(Np), np.zeros(Np)
u,ud,udd = u0,ud0,udd0

u_last = ud0_last*T

ff = -Poids_x

fp_first = np.zeros(Np)
fp = np.zeros(Nt)
p = 0
while p<Nt :
    if T[p] >= T_Fp :
        fp[p] = P_gli[0]*np.sin(np.pi*(T[p]-T_Fp))
    p+=1

# paramètres d'affichage du problème
U_plot,Ud_plot,Udd_plot = [u0],[ud0],[udd0]
ff_plot, f_plot, f_trial_plot = [], [], []
Ff_plot, F_plot, F_trial_plot = [], [], []

Niter_implicite,ErrResi = [],[]

# matrix and coefficient
alpha_HHT = -1/3

mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

Adh = [False for k in range(Np)]
Adh_first = [True for k in range(Np)]
Adh_last = Adh.copy()
Adh_list = []
Sign_list = []
ind_zoom = 0
sign_last = np.sign(ud0)
f_trial_n = np.zeros(Np)

# first step for the external load
# connection with the wall moving at ud0_last speed
k_last = K[-1][-1]
fe_last_list = [0 for k in range(Np-1)]
fe_last_list.append(-k_last*(u[-1] - u_last[0]))
fe_last = np.array(fe_last_list)
del fe_last_list

# elastic force due to the stiffness between blocs
fe = -np.dot(K,u)

# external loading on the first blocs
fp_first_list = []
fp_first_list.append(fp[0])
for k in range(Np-1):
    fp_first_list.append(0)
fp_first = np.array(fp_first_list)
del fp_first_list

p_n = fe + fe_last + Poids_x + fp_first + ff

for n in range(Nt-1):
    
    Adh = Adh_last.copy()
    
    un = u
    udn = ud
    uddn = udd
    
    u_i = u + dt*ud + 0.5*(dt**2)*udd
    ud_i = ud + dt*udd
    udd_i = udd
    
    ms_i = M + mk*K
    
    ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
    fe_last_list = [0 for k in range(Np-1)]
    fe_last_list.append(-k_last*(u_i[-1] - u_last[n+1]))
    fe_last = np.array(fe_last_list)
    del fe_last_list
    
    fe = -np.dot(K,np.transpose(u_i))
    
    fp_first_list = []
    fp_first_list.append(fp[n+1])
    for k in range(Np-1):
        fp_first_list.append(0)
    fp_first = np.array(fp_first_list)
    del fp_first_list
    
    f_trial_i = fe + fe_last + Poids_x + fp_first
    
    ## définition des états de glissement et adhérence au premier bloc de Coulomb
    for i in range(Np):
        if Adh_last[i] == False :
            if ud_i[i]==0 or udn[i]*ud_i[i]<0 or np.abs(f_trial_i[i]) <= P_gli[i]:
                Adh[i] = True
        elif Adh_last[i] == True :
            if np.abs(f_trial_i[i]) > P_gli[i]:
                Adh[i] = False
                if ud_i[i] != 0:
                    sign_last[i] = np.sign(ud_i[i])
            if f_trial_i[i]*f_trial_n[i] < 0:
                sign_last[i] = np.sign(f_trial_i[i])
    
    
    ## calcul des forces Coulomb et calcul des résultantes
    ff_alpha_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
    p_alpha_t_list = []
    
    for i in range(Np):
        if Adh[i] == True : 
            ff_alpha_list.append(-np.abs(f_trial_i[i])*sign_last[i])
            p_alpha_t_list.append(0)
        if Adh[i] == False : 
            ff_alpha_list.append(-P_gli[i]*sign_last[i])
            p_alpha_t_list.append(ff_alpha_list[i] + f_trial_i[i])
    
    ff_alpha = np.array(ff_alpha_list)
    p_alpha_t = np.array(p_alpha_t_list)
    del ff_alpha_list,p_alpha_t_list
    
    # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
    r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(M,udd_i)
    delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
    u_i1 = u_i + uiai*delta_uddn_i
    ud_i1 = ud_i + udiai*delta_uddn_i
    udd_i1 = udd_i + delta_uddn_i
    
    for i in range(Np):
        if Adh[i] == True :
            u_i1[i] = u_i[i]
            ud_i1[i] = 0
            udd_i1[i] = 0
    
    iter = 0
    
    while iter<10 and np.max(np.abs(udd_i1-udd_i))>eps_HHT :
        
        Adh_last = Adh.copy()
        
        ms_i = M + mk*K
        
        ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
        fe_last_list = [0 for i in range(Np-1)]
        fe_last_list.append(-K[1][1]*(u_i1[1] - u_last[n+1]))
        fe_last = np.array(fe_last_list)
        del fe_last_list
        
        fe = -np.dot(K,np.transpose(u_i1))
        
        fp_first_list = []
        fp_first_list.append(fp[n+1])
        for k in range(Np-1):
            fp_first_list.append(0)
        fp_first = np.array(fp_first_list)
        del fp_first_list
        
        f_trial_i = fe + fe_last + Poids_x + fp_first
        
        ## définition des états de glissement et adhérence au premier bloc de Coulomb
        for i in range(Np):
            if Adh_last[i] == False :
                if ud_i1[i]==0 or udn[i]*ud_i1[i]<0 or np.abs(f_trial_i[i]) <= P_gli[i]:
                    Adh[i] = True
            elif Adh_last[i] == True :
                if np.abs(f_trial_i[i]) > P_gli[i]:
                    Adh[i] = False
                    if ud_i1[i] != 0:
                        sign_last[i] = np.sign(ud_i1[i])
                if f_trial_i[i]*f_trial_n[i] < 0:
                    sign_last[i] = np.sign(f_trial_i[i])
        
        
        ## calcul des forces Coulomb et calcul des résultantes
        ff_alpha_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
        p_alpha_t_list = []
        
        for i in range(Np):
            if Adh[i] == True : 
                ff_alpha_list.append(-np.abs(f_trial_i[i])*sign_last[i])
                p_alpha_t_list.append(0)
            if Adh[i] == False : 
                ff_alpha_list.append(-P_gli[i]*sign_last[i])
                p_alpha_t_list.append(ff_alpha_list[i] + f_trial_i[i])
        
        ff_alpha = np.array(ff_alpha_list)
        p_alpha_t = np.array(p_alpha_t_list)
        del ff_alpha_list,p_alpha_t_list
        
        u_i = u_i1
        ud_i = ud_i1
        udd_i = udd_i1
        
        # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
        r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(M,udd_i)
        delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
        u_i1 = u_i + uiai*delta_uddn_i
        ud_i1 = ud_i + udiai*delta_uddn_i
        udd_i1 = udd_i + delta_uddn_i
        
        for i in range(Np):
            if Adh[i] == True :
                u_i1[i] = u_i[i]
                ud_i1[i] = 0
                udd_i1[i] = 0
    
        iter +=1
    
    if iter==10 :
        print('Alerte non convergence au pas ' + str(n))
    
    ff = ff_alpha.copy()
    p_n = p_alpha_t.copy()
    f_trial_n = f_trial_i.copy()
    
    Adh_last = Adh.copy()
    del Adh
    
    en = np.abs(ud_i1 - ud_i) #residual error
    
    # state vector at n+1
    u = u_i1
    ud = ud_i1
    udd = udd_i1
    
    if n%1==0 :
        
        U_plot.append( u )
        Ud_plot.append( ud )
        Udd_plot.append( udd )
        
        Ff_plot.append( ff )
        F_plot.append( p_alpha_t )
        F_trial_plot.append( f_trial_i )
        Adh_list.append( Adh_last )
        Sign_list.append( sign_last ) 
        
        Niter_implicite.append(iter)
        ErrResi.append(en)

save_figure_path = 'C:\\users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_11_13_Resultats_2blocs_Coulomb'
os.chdir(save_figure_path)

import matplotlib
matplotlib.rcParams['figure.figsize'] = 35,18
matplotlib.rcParams.update({'font.size': 20})
title_font = {'size':'30'}
label_font = {'size':'30'}

T_plot = np.linspace(0,T_end,len(U_plot))

Cc_label = 'Fp_dt000314'
plt.figure(1)
plt.plot(T_plot,np.array(U_plot)*1e+3,linewidth=2)
# plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([u_last[0],u_last[-1]])*1e+3,'b--',linewidth=2)
# # plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
# plt.xlabel('Temps ($s$)',label_font)
# plt.ylabel('Déplacement ($mm$)',label_font)
# plt.title("Déplacement du modèle 2 blocs avec $\mu = tan( alpha ) $",label_font)
# plt.grid(True)
# # plt.axis([0,20,-5,5])
# plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0 .T$'],loc='best')
# # plt.savefig('DeplacementTestCc' + Cc_label + '.png')
# plt.savefig('Deplacement_' + Cc_label + '.pdf')

plt.figure(2)
plt.plot(T_plot,np.array(Ud_plot)*1e+3,linewidth=2)
# plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*ud0_last*1e+3,'y--',linewidth=2)
# # plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
# plt.xlabel('Temps ($s$)',label_font)
# plt.ylabel('Vitesse ($mm/s$)',label_font)
# plt.title("Vitesse du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
# plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0$'],loc='best')
# plt.grid(True)
# # plt.axis([0,20,-5,5])
# # plt.savefig('VitesseTestCc' + Cc_label + '.png')
# plt.savefig('Vitesse_' + Cc_label + '.pdf')

fp_plot_list = []
for k in range(len(T_plot)):
    if k < len(T_plot) :
        fp_plot_list.append(fp[k])
    else :
        fp_plot_list.append(fp[-1])
fp_plot = np.array(fp_plot_list)

plt.figure(3)
plt.plot(T_plot[1:],np.array(Ff_plot)*1e-6,linewidth=2)
# plt.plot(T_plot,fp_plot*1e-6,'y--',linewidth=2)
# # plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
# plt.plot(np.array([T_plot[0],T_plot[-1]]),-np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
# plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
# plt.xlabel('Temps ($s$)',label_font)
# plt.ylabel('Force de frottement ($MN/m$)',label_font)
# plt.title("Frottement du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
# plt.legend(['Bloc Coulomb','Bloc Weertman','Force de perturbation','Limites $-P_{gli}$ et $P_{gli}$ de glissement'],loc='best')
# plt.grid(True)
# # plt.savefig('FrottementTestCc' + Cc_label + '.png')
# plt.savefig('Frottement_' + Cc_label + '.pdf')

plt.figure(4)
plt.plot(T_plot[1:],np.array(F_plot)*1e-6,linewidth=2)
# # plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
# plt.xlabel('Temps ($s$)',label_font)
# plt.ylabel('Resultante ($MN/m$)',label_font)
# plt.title("Résultante du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
# plt.grid(True)
# # plt.axis([0,20,-10,10])
# plt.legend(['Bloc Coulomb','Bloc Weertman'],loc='best')
# # plt.savefig('ResultanteTestCc' + Cc_label + '.png')
# plt.savefig('Resultante_' + Cc_label + '.pdf')

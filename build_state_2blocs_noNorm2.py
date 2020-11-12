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
Ltot = 200/np.cos(alpha_ground)
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Xt = -1300
Xc = -1400

# définition de la géométrie du glacier
X_down = np.array([Xt,Xc,Xt-Ltot*np.cos(alpha_ground)])
dl_l = np.array([X_down[0]-X_down[1],X_down[1]-X_down[2]])
dl_min = np.min(dl_l)

M_bloc = rho_glace*H*Ltot/2
Poids_x = -np.sin(alpha_ground)*M_bloc*g* np.array([1,1])
Poids_y = -np.cos(alpha_ground)*M_bloc*g* np.array([1,1])
P_gli = Cc*np.abs(Poids_y[0])

# Steady sliding
Cw = (np.abs(Poids_x[1]))/(dl_l[1]*(np.abs(ud0_last))**(m))
ud_reg = np.abs(ud0_last)
Cl = Cw*(ud_reg)**(m-1)

eps_HHT = P_gli/(100*M_bloc)

K = E*H/(Ltot/2)*np.array([[1,-1],[-1,1]])
M = M_bloc*np.array([[1,0],[0,1]])
T_end = 5
T_Fp = 2

# discretisation des variables
dt = 0.1*dl_min*np.sqrt(rho_glace/E)
Nt = int(T_end/dt) + 1
T = np.linspace(0,T_end,Nt)


# definition des variables du problème
u0,ud0,udd0 = np.array([0,0]), np.array([ud0_last,ud0_last]), np.array([0,0])
u,ud,udd = u0,ud0,udd0

u_last = ud0_last*T

ff = -Poids_x
fe_last = np.zeros(2)

fp_first = np.zeros(2)
fp = np.zeros(Nt)
p = 0
while p<Nt :
    if T[p] >= T_Fp :
        fp[p] = P_gli*np.sin(np.pi*(T[p]-T_Fp))
    p+=1

# paramètres d'affichage du problème
U_plot,Ud_plot,Udd_plot = [u0],[ud0],[udd0]
ff_plot, f_plot, f_trial_plot = [], [], []
Ff_plot, F_plot, F_trial_plot = [], [], []

Niter_implicite,ErrResi = [],[]

# matrix and coefficient
alpha_HHT = -1/3

mc = (1+alpha_HHT)*(0.5-alpha_HHT)*dt
mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

Adh = [False,False]
Adh_first = [True,True]
Adh_last = Adh.copy()
Adh_list = []
Sign_list = []
ind_zoom = 0
sign_last = np.sign(ud0)

fe_last[1] = -K[1][1]*(u[1] - u_last[0])
fe = -np.dot(K,u)

fp_first[0] = fp[0]

p_n = fe + fe_last + Poids_x + fp_first + ff

for n in range(Nt-1):
    
    Adh = Adh_last.copy()
    
    un = u
    udn = ud
    uddn = udd
    
    u_i = u + dt*ud + 0.5*(dt**2)*udd
    ud_i = ud + dt*udd
    udd_i = udd
    
    C_i = np.zeros((2,2))
    if np.abs(ud_i[1]) > ud_reg :
        C_i[1][1] = Cw*m*dl_l[1]*((np.abs(ud_i[1]))**(m-1))
    else :
        C_i[1][1] = Cl*dl_l[1]
    
    ms_i = M + mc*C_i + mk*K
    
    ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
    fe_last[1] = -K[1][1]*(u_i[1] - u_last[n+1])
    fe = -np.dot(K,np.transpose(u_i))
    fp_first[0] = fp[n+1]
    
    f_trial_i = fe + fe_last + Poids_x + fp_first
    
    ## définition des états de glissement et adhérence au premier bloc de Coulomb
    if Adh_last[0] == False :
        if ud_i[0]==0 or udn[0]*ud_i[0]<0 or np.abs(f_trial_i[0]) <= P_gli: # passage de glissement à adhérence
            print('Test adhérence completed au pas :' +str(n) + ' à la première itération')
            Adh[0] = True
    elif Adh_last[0] == True :
        if np.abs(f_trial_i[0]) > P_gli : # passage d'adhérence à glissement
            Adh[0] = False
            print('Test glissement completed au pas :' +str(n) + ' à la première itération')
            if ud_i[0] != 0:
                sign_last[0] = np.sign(ud_i[0])
            else :
                sign_last[0] = np.sign(f_trial_i[0])
    
    ## frottement de Weertman pour le deuxième bloc
    ff_weert_list = [0]
    if np.abs(ud_i[1])>ud_reg :
        ff_weert_list.append(-np.sign(ud_i[1])*Cw*dl_l[1]*((np.abs(ud_i[1]))**m))
    else :
        ff_weert_list.append(-Cl*dl_l[1]*ud_i[1])
    ff_weert = np.array(ff_weert_list)
    
    ## calcul des forces Coulomb et calcul des résultantes
    ff_coul_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
    p_alpha_t_list = []
    
    if Adh[0] == True : 
        ff_coul_list.append(f_trial_i[0]*sign_last[0])
        p_alpha_t_list.append(0)
    if Adh[0] == False : 
        ff_coul_list.append(-P_gli*sign_last[0])
        p_alpha_t_list.append(ff_coul_list[0] + f_trial_i[0])
    
    ff_coul_list.append(0)
    ff_coul = np.array(ff_coul_list)
    
    ff_alpha = ff_coul + ff_weert
    
    p_alpha_t_list.append(f_trial_i[1] + ff_weert[1])
    p_alpha_t = np.array(p_alpha_t_list)
    
    p_alpha_t[1] = f_trial_i[1] + ff_alpha[1]
    
    # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
    r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(M,udd_i)
    delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
    u_i1 = u_i + uiai*delta_uddn_i
    ud_i1 = ud_i + udiai*delta_uddn_i
    udd_i1 = udd_i + delta_uddn_i
    
    if Adh[0] == True :
        u_i1[0] = u_i[0]
        ud_i1[0] = 0
        udd_i1[0] = 0
    
    i = 0
    
    while i<10 and np.max(np.abs(udd_i1-udd_i))>eps_HHT :
        
        Adh_last = Adh.copy()
        
        c_i = np.zeros((2,2))
        if np.abs(ud_i[1]) > ud_reg :
            c_i[1][1] = Cw*m*dl_l[1]*((np.abs(ud_i[1]))**(m-1))
        else :
            c_i[1][1] = Cl*dl_l[1]
        
        ms_i = M + mc*c_i + mk*K
        
        ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
        fe_last[1] = -K[1][1]*(u_i1[1] - u_last[n+1])
        fe = -np.dot(K,np.transpose(u_i1))
        fp_first[0] = fp[n+1]
        
        f_trial_i = fe + fe_last + Poids_x + fp_first
        
        ## définition des états de glissement et adhérence au premier bloc de Coulomb
        if Adh[0] == False :
            if ud_i1[0]==0 or udn[0]*ud_i1[0]<0 or np.abs(f_trial_i[0]) <= P_gli: # passage de glissement à adhérence
                print('Test adhérence completed au pas :' +str(n) + ' itération '+ str(i+1))
                Adh[0] = True
        elif Adh[0] == True :
            if np.abs(f_trial_i[0]) > P_gli : # passage d'adhérence à glissement
                print('Test glissement completed au pas :' +str(n) + ' itération '+ str(i+1))
                Adh[0] = False
                if ud_i1[0] != 0:
                    sign_last[0] = np.sign(ud_i1[0])
                else :
                    sign_last[0] = np.sign(f_trial_i[0])
        
        ## frottement de Weertman pour le deuxième bloc
        ff_weert_list = [0]
        if np.abs(ud_i1[1])>ud_reg :
            ff_weert_list.append(-np.sign(ud_i1[1])*Cw*dl_l[1]*((np.abs(ud_i1[1]))**m))
        else :
            ff_weert_list.append(-Cl*dl_l[1]*ud_i1[1])
        ff_weert = np.array(ff_weert_list)
        
        ## calcul des forces Coulomb et calcul des résultantes
        ff_alpha_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
        p_alpha_t = np.zeros(2)
        
        if Adh[0] == True : 
            ff_alpha_list.append(f_trial_i[0]*sign_last[0])
            p_alpha_t[0] = 0
        if Adh[0] == False : 
            ff_alpha_list.append(-P_gli*sign_last[0])
            p_alpha_t[0] = ff_alpha[0] + f_trial_i[0]
        
        ff_alpha_list.append(ff_weert[1])
        ff_alpha = np.array(ff_alpha_list)
        
        p_alpha_t[1] = f_trial_i[1] + ff_alpha[1]
        
        u_i = u_i1
        ud_i = ud_i1
        udd_i = udd_i1
        
        # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
        r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(M,udd_i)
        delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
        u_i1 = u_i + uiai*delta_uddn_i
        ud_i1 = ud_i + udiai*delta_uddn_i
        udd_i1 = udd_i + delta_uddn_i
        
        if Adh[0] == True :
            u_i1[0] = u_i[0]
            ud_i1[0] = 0
            udd_i1[0] = 0
    
        i +=1
    
    if i==10 :
        print('Alerte non convergence au pas ' + str(n))
    
    ff = ff_alpha.copy()
    p_n = p_alpha_t.copy()
    
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
        
        Niter_implicite.append(i)
        ErrResi.append(en)

save_figure_path = 'C:\\users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_10_19_CoulWeert_2blocs_test_CcTan'
os.chdir(save_figure_path)

import matplotlib
matplotlib.rcParams['figure.figsize'] = 35,18
matplotlib.rcParams.update({'font.size': 20})
title_font = {'size':'30'}
label_font = {'size':'30'}

T_plot = np.linspace(0,T_end,len(U_plot))

Cc_label = 'Tan_init3_dt00314_noNorm'
plt.figure(1)
plt.plot(T_plot,np.array(U_plot)*1e+3,linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([u_last[0],u_last[-1]])*1e+3,'b--',linewidth=2)
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
plt.plot(T_plot,np.array(Ud_plot)*1e+3,linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*ud0_last*1e+3,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Vitesse ($mm/s$)',label_font)
plt.title("Vitesse du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0$'],loc='best')
plt.grid(True)
# plt.axis([0,20,-5,5])
plt.savefig('VitesseTestCc' + Cc_label + '.png')
plt.savefig('VitesseTestCc' + Cc_label + '.pdf')

fp_plot_list = []
for k in range(len(T_plot)):
    if k < len(T_plot) :
        fp_plot_list.append(fp[k])
    else :
        fp_plot_list.append(fp[-1])
fp_plot = np.array(fp_plot_list)

plt.figure(3)
plt.plot(T_plot[1:],np.array(Ff_plot)*1e-6,linewidth=2)
plt.plot(T_plot,fp_plot*1e-6,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),-np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Force de frottement ($MN/m$)',label_font)
plt.title("Frottement du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.legend(['Bloc Coulomb','Bloc Weertman','Force de perturbation','Limites $-P_{gli}$ et $P_{gli}$ de glissement'],loc='best')
plt.grid(True)
plt.savefig('FrottementTestCc' + Cc_label + '.png')
plt.savefig('FrottementTestCc' + Cc_label + '.pdf')

plt.figure(4)
plt.plot(T_plot[1:],np.array(F_plot)*1e-6,linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Resultante ($MN/m$)',label_font)
plt.title("Résultante du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.grid(True)
# plt.axis([0,20,-10,10])
plt.legend(['Bloc Coulomb','Bloc Weertman'],loc='best')
plt.savefig('ResultanteTestCc' + Cc_label + '.png')
plt.savefig('ResultanteTestCc' + Cc_label + '.pdf')

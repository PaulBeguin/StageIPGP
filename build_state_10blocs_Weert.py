
# importation des modules
import numpy as np
import matplotlib.pyplot as plt
import os


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
m = 1
# Cc = np.tan(alpha_ground)
Cc = np.tan(alpha_ground)
type_law = 1 #0=Weertman, 1=Coulomb

# Geometric
H = 800
Ltot = 200/np.cos(alpha_ground)
h_im = H * (rho_glace/rho_eau) * np.cos(alpha_ground)
Xt = -1300
Xc = -1400
Np = 3

## geometric definition
dl_middle = Ltot/(Np-1)

# X coordinates
X_down_list = []
X_down_list.append(Xt)
X_down_list.append(Xt - dl_middle/2)
for i in range(Np-2):
    X_down_list.append(X_down_list[-1]-dl_middle)
X_down_list.append(X_down_list[-1]-dl_middle/2)
X_down = np.array(X_down_list)

# length of the blocs
dl_l_list = [X_down_list[i]-X_down_list[i+1] for i in range(Np)]
dl_l = np.array(dl_l_list)

# lenght of the elastic bonds
dl_liai_list = []
dl_liai_list.append(dl_l_list[0] + dl_l_list[1]/2)
for i in range(Np-3):
    dl_liai_list.append((dl_l_list[i+1] + dl_l_list[i+2])/2)
dl_liai_list.append(dl_l_list[Np-2]/2 + dl_l_list[Np-1])
dl_liai = np.array(dl_liai_list)

dl_min = np.min(dl_l)

# stiffness matrix
K = np.zeros((Np,Np))
ke = E*H*np.array([[1,-1],[-1,1]])
for i in range(Np-1):
    K[i:i+2,i:i+2] = K[i:i+2,i:i+2] + ke/dl_liai[i]

# mass matrix
M_bloc_list = []
for i in range(Np):
    M_bloc_list.append(rho_glace*H*dl_l[i])
M_bloc = np.array(M_bloc_list)
M = np.diag(M_bloc)

# weight array
Poids_x = -np.sin(alpha_ground)*g*M_bloc
Poids_y = -np.cos(alpha_ground)*g*M_bloc

# steady sliding
Cw = (np.abs(Poids_x[1]))/(dl_l[1]*(np.abs(ud0_last))**(m))
ud_reg = np.abs(ud0_last)
Cl = Cw*(ud_reg)**(m-1)

# time parameter 
T_end = 10
T_Fp = 2

dt = 0.1*dl_min*np.sqrt(rho_glace/E)
Nt = int(T_end/dt) + 1
T = np.linspace(0,T_end,Nt)

# definition des variables du problème
u0,ud0,udd0 = np.zeros(Np), ud0_last*np.ones(Np), np.zeros(Np)
u,ud,udd = u0,ud0,udd0

u_last = ud0_last*T

ff = -Poids_x

fp = np.zeros(Nt)
p = 0
while p<Nt :
    if T[p] >= T_Fp :
        fp[p] = Poids_x[0]*np.sin(np.pi*(T[p]-T_Fp))
    p+=1

# paramètres d'affichage du problème
U_plot,Ud_plot,Udd_plot = [u0],[ud0],[udd0]
Ff_plot, F_plot, F_sis = [], [], []

Niter_implicite,ErrResi = [],[]

# matrix and coefficient
alpha_HHT = -1/3
eps_HHT = np.min(np.abs(Poids_x))/(100*np.min(M_bloc))

mc = (1+alpha_HHT)*(0.5-alpha_HHT)*dt
mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

## first step
# elastic force at the last bloc
fe_last_list = []
for i in range(Np-1):
    fe_last_list.append(0)
fe_last_list.append( -K[-1][-1]*(u[-1] - u_last[0]))
fe_last = np.array(fe_last_list)

# elastic force
fe = -np.dot(K,u)

# loading on the first bloc
fp_first_list = []
fp_first_list.append(fp[0])
for i in range(Np-1):
    fp_first_list.append(0)
fp_first = np.array(fp_first_list)

p_n = fe + fe_last + Poids_x + fp_first + ff

## start of time integration
for n in range(Nt-1):
    
    un = u
    udn = ud
    uddn = udd
    
    u_i = u + dt*ud + 0.5*(dt**2)*udd
    ud_i = ud + dt*udd
    udd_i = udd
    
    u_i1,ud_i1,udd_i1 = np.zeros(Np),np.zeros(Np),np.zeros(Np)
    
    # exernal load on the first bloc at n+1 step of time
    fp_first_list = []
    fp_first_list.append(fp[n+1])
    for k in range(Np-1):
        fp_first_list.append(0)
    fp_first = np.array(fp_first_list)
    
    iter = 0
    
    while iter==0 or (iter < 10 and np.max(udd_i1 - udd_i)>eps_HHT ) :
        
        if iter > 0 :
            u_i = u_i1
            ud_i = ud_i1
            udd_i = udd_i1
            
        # damping matrix
        C_i_list = []
        for k in range(Np):
            if np.abs(ud_i[k]) >= ud_reg :
                C_i_list.append( Cw*m*dl_l[k]*((np.abs(ud_i[k]))**(m-1)) )
            else :
                C_i_list.append( Cl*dl_l[k] )
        C_i = np.diag(C_i_list)
        
        ms_i = M + mc*C_i + mk*K
        
        ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
        fe_last_list = []
        for k in range(Np-1):
            fe_last_list.append(0)
        fe_last_list.append(-K[-1][-1]*(u_i[-1] - u_last[n+1]))
        fe_last = np.array(fe_last_list)
        
        # elastic bonds
        fe = -np.dot(K,np.transpose(u_i))
        
        # Weertman's sliding on all blocs
        ff_weert_list = []
        for k in range(Np):
            if np.abs(ud_i[k])>=ud_reg :
                ff_weert_list.append(-np.sign(ud_i[k])*Cw*dl_l[k]*((np.abs(ud_i[k]))**m))
            else :
                ff_weert_list.append(-Cl*dl_l[k]*ud_i[k])
        ff_weert = np.array(ff_weert_list)
        
        p_alpha_t = fe + fe_last + Poids_x + fp_first + ff_weert
        
        # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
        r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(M,udd_i)
        delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
        
        u_i1 = u_i + uiai*delta_uddn_i
        ud_i1 = ud_i + udiai*delta_uddn_i
        udd_i1 = udd_i + delta_uddn_i
        
        iter += 1
    
    if iter == 10 :
        print('Alerte non convergence au pas ' + str(n))
    
    ff = ff_weert.copy()
    p_n = p_alpha_t.copy()
    
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
        F_plot.append( p_alpha_t - fp_first )
        F_sis.append ( np.sum(ff) + np.sum(Poids_x) + fe_last[-1] )
        
        Niter_implicite.append(iter)
        ErrResi.append(en)

save_figure_path = 'C:\\users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_11_16_Weert_10blocs'
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
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([u_last[0],u_last[-1]])*1e+3,'b--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Déplacement ($mm$)',label_font)
plt.title("Déplacement du modèle 2 blocs avec $\mu = tan( alpha ) $",label_font)
plt.grid(True)
# plt.axis([0,20,-5,5])
# plt.savefig('DeplacementTestCc' + Cc_label + '.png')
plt.savefig('Deplacement_' + Cc_label + '.pdf')

plt.figure(2)
plt.plot(T_plot,np.array(Ud_plot)*1e+3,linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*ud0_last*1e+3,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Vitesse ($mm/s$)',label_font)
plt.title("Vitesse du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.grid(True)
# plt.axis([0,20,-5,5])
# plt.savefig('VitesseTestCc' + Cc_label + '.png')
plt.savefig('Vitesse_' + Cc_label + '.pdf')

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
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Force de frottement ($MN/m$)',label_font)
plt.title("Frottement du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.grid(True)
# plt.savefig('FrottementTestCc' + Cc_label + '.png')
plt.savefig('Frottement_' + Cc_label + '.pdf')

plt.figure(4)
plt.plot(T_plot[1:],np.array(F_plot)*1e-6,linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Resultante ($MN/m$)',label_font)
plt.title("Résultante du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.grid(True)
# plt.axis([0,20,-10,10])
# plt.savefig('ResultanteTestCc' + Cc_label + '.png')
plt.savefig('Resultante_' + Cc_label + '.pdf')


# sismic force
plt.figure(5)
plt.plot(T_plot[1:],F_sis)
plt.plot(T_plot,-fp)

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

# Steady sliding
Cw = (rho_glace*H*g*np.sin(alpha_ground))/((np.abs(Ud0_last))**(m))
Ud_reg = 0.1*np.abs(Ud0_last)
Cl = Cw*(Ud_reg)**(m-1)

# définition de la géométrie du glacier
X_down = np.array([Xt,Xc,Xt-Ltot*np.cos(alpha_ground)])
dl_l = np.array([X_down[0]-X_down[1],X_down[1]-X_down[2]])
dl_min = np.min(dl_l)

M_bloc = rho_glace*H*Ltot/2
Poids_x = -np.sin(alpha_ground)*M_bloc*g* np.array([1,1])
Poids_y = -np.cos(alpha_ground)*M_bloc*g* np.array([1,1])
P_gli = Cc*np.abs(Poids_y[0])

eps_HHT = P_gli/(100*M_bloc)

K = E*H/(Ltot/2)*np.array([[1,-1],[-1,1]])
M = M_bloc*np.array([[1,0],[0,1]])
<<<<<<< HEAD
T_end = 5
T_Fp = 2
=======
T_end = 17
T_Fp = 15
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc

# normalisation du problème
L_norm = dl_min
M_norm = rho_glace*(dl_min**3)
F_norm = P_gli
T_norm = dl_min*dl_min*np.sqrt(rho_glace/P_gli)
Cw0 = F_norm/(L_norm*(L_norm/T_norm)**(m))
Cl0 = F_norm/(L_norm*(L_norm/T_norm))

# discretisation des variables
<<<<<<< HEAD
dt = 0.0001*dl_min*np.sqrt(rho_glace/E)
Nt = int(T_end/dt) + 1

# T_plot = np.linspace(0,T_end,Nt)
=======
dt = 0.001*dl_min*np.sqrt(rho_glace/E)
Nt = int(T_end/dt) + 1

T_plot = np.linspace(0,T_end,Nt)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
t_norm = np.linspace(0,T_end/T_norm,Nt)

x_down = X_down/L_norm
dl_norm = dl_l/L_norm

ud0_last = Ud0_last/(L_norm/T_norm)
ud_reg = Ud_reg/(L_norm/T_norm)

m_norm = M/M_norm
<<<<<<< HEAD
# k_norm = K/(E*Ltot**3)
k_norm = K/(F_norm/L_norm)
=======
k_norm = K/(E*Ltot**3)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc

cw_norm = Cw/Cw0
cl_norm = Cl/Cl0

p_gli = P_gli/F_norm
poids_x = Poids_x/F_norm
poids_y = Poids_y/F_norm

# definition des variables du problème
<<<<<<< HEAD
u0,ud0,udd0 = np.array([0,0]), np.array([-1,-1]), np.array([0,0])
=======
u0,ud0,udd0 = np.array([0,0]), np.array([0,0]), np.array([0,0])
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
u,ud,udd = u0,ud0,udd0

u_last = ud0_last*t_norm

<<<<<<< HEAD
ff = np.ones(2)
=======
ff = np.zeros(2)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
fe_last = np.zeros(2)

fp_first = np.zeros(2)
fp = np.zeros(Nt)
t_Fp = T_Fp/T_norm
p = 0
while p<Nt :
    if t_norm[p] >= t_Fp :
        fp[p] = p_gli*np.sin(np.pi*(t_norm[p]-t_Fp)*T_norm)
    p+=1

# paramètres d'affichage du problème
u_plot,ud_plot,udd_plot = [u0],[ud0],[udd0]
U_plot,Ud_plot,Udd_plot = [u0*Ltot],[ud0*(Ltot/T_norm)],[udd0*(Ltot/(T_norm*T_norm))]
ff_plot, f_plot, f_trial_plot = [], [], []
Ff_plot, F_plot, F_trial_plot = [], [], []

Niter_implicite,ErrResi = [],[]

# matrix and coefficient
<<<<<<< HEAD
alpha_HHT = -1/3
=======
alpha_HHT = -0.005
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc

mc = (1+alpha_HHT)*(0.5-alpha_HHT)*dt
mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
uiai = 0.25*((1-alpha_HHT)*dt)**2
udiai = (0.5-alpha_HHT)*dt
p_alpha = -alpha_HHT
p_alpha1 = (1+alpha_HHT)

Adh = [False,False]
Adh_first = [True,True]
<<<<<<< HEAD
Adh_last = Adh.copy()
Adh_list = []
Sign_list = []
=======
Adh_list = []
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
ind_zoom = 0
sign_last = np.array([1,1])

fe_last[1] = -k_norm[1][1]*(u[1] - u_last[0])
fe = -np.dot(k_norm,u)

fp_first[0] = fp[0]

p_n = fe + fe_last + poids_x + fp_first + ff

for n in range(Nt-1):
    
<<<<<<< HEAD
    Adh = Adh_last.copy()
    
=======
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    un = u
    udn = ud
    uddn = udd
    
    u_i = u + dt*ud + 0.5*(dt**2)*udd
    ud_i = ud + dt*udd
    udd_i = udd
    
    c_i = np.zeros((2,2))
    if np.abs(ud_i[1]) > ud_reg :
<<<<<<< HEAD
        c_i[1][1] = (Cw*m*dl_l[1]*(((Ltot/T_norm)*np.abs(ud_i[1]))**(m-1)))/M_norm
    else :
        c_i[1][1] = (Cl*dl_l[1])/M_norm
=======
        c_i[1][1] = cw_norm*m*dl_norm[1]*(np.abs(ud_i[1])**(m-1))
    else :
        c_i[1][1] = cl_norm*dl_norm[1]
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    
    ms_i = m_norm + mc*c_i + mk*k_norm
    
    ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
    fe_last[1] = -k_norm[1][1]*(u_i[1] - u_last[n+1])
    fe = -np.dot(k_norm,np.transpose(u_i))
    fp_first[0] = fp[n+1]
    
    f_trial_i = fe + fe_last + poids_x + fp_first
    
    ## définition des états de glissement et adhérence au premier bloc de Coulomb
<<<<<<< HEAD
    if Adh_last[0] == False :
        if ud_i[0]==0 or udn[0]*ud_i[0]<0 or np.abs(f_trial_i[0]) <= p_gli: # passage de glissement à adhérence
            print('Test adhérence completed au pas :' +str(n) + ' à la première itération')
            Adh[0] = True
    elif Adh_last[0] == True :
        if np.abs(f_trial_i[0]) > p_gli : # passage d'adhérence à glissement
            Adh[0] = False
            print('Test glissement completed au pas :' +str(n) + ' à la première itération')
            if ud_i[0] != 0:
=======
    if Adh[0] == False :
        if ud_i[0]==0.0 or udn[0]*ud_i[0]<0 or np.abs(f_trial_i[0]) <= p_gli: # passage de glissement à adhérence
            Adh[0] = True
    elif Adh[0] == True :
        if np.abs(f_trial_i[0]) > p_gli : # passage d'adhérence à glissement
            Adh[0] = False
            if ud_i[0] != 0.0:
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
                sign_last[0] = np.sign(ud_i[0])
            else :
                sign_last[0] = np.sign(f_trial_i[0])
    
    ## frottement de Weertman pour le deuxième bloc
<<<<<<< HEAD
    ff_weert_list = [0]
    if np.abs(ud_i[1])>ud_reg :
        ff_weert_list.append(-(np.sign(ud_i[1])*Cw*dl_l[1]*(((L_norm/T_norm)*np.abs(ud_i[1]))**(m)))/F_norm)
    else :
        ff_weert_list.append(-(Cl*dl_l[1]*(L_norm/T_norm)*ud_i[1])/F_norm)
    ff_weert = np.array(ff_weert_list)
    
    ## calcul des forces Coulomb et calcul des résultantes
    ff_alpha_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
    p_alpha_t = np.zeros(2)
    
    if Adh[0] == True : 
        ff_alpha_list.append(f_trial_i[0]*sign_last[0])
        p_alpha_t[0] = 0
    if Adh[0] == False : 
        ff_alpha_list.append(-p_gli*sign_last[0])
        p_alpha_t[0] = ff_alpha_list + f_trial_i[0]
    
    ff_alpha_list.append(ff_weert_list[1])
    ff_alpha = np.array(ff_alpha_list)
=======
    ff_weert = np.array([0,0])
    if np.abs(ud_i[1])>ud_reg :
        ff_weert[1]=-np.sign(ud_i[1])*cw_norm*dl_norm[1]*(np.abs(ud_i[1])**(m))
    else :
        ff_weert[1]=-cl_norm*dl_norm[1]*ud_i[1]
    
    ## calcul des forces Coulomb et calcul des résultantes
    ff_alpha = ff_weert # initialisation avec le frottement aux blocs de Weertman seulement
    p_alpha_t = np.zeros(2)
    
    if Adh[0] == True : 
        ff_alpha[0] = -f_trial_i[0]*sign_last[0]
        p_alpha_t[0] = 0
    if Adh[0] == False : 
        ff_alpha[0] = -p_gli*sign_last[0]
        p_alpha_t[0] = ff_alpha[0] + f_trial_i[0]
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    
    p_alpha_t[1] = f_trial_i[1] + ff_alpha[1]
    
    # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
    r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(m_norm,udd_i)
    delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
    u_i1 = u_i + uiai*delta_uddn_i
    ud_i1 = ud_i + udiai*delta_uddn_i
    udd_i1 = udd_i + delta_uddn_i
    
    if Adh[0] == True :
        u_i1[0] = u_i[0]
<<<<<<< HEAD
        ud_i1[0] = 0
        udd_i1[0] = 0
=======
        ud_i1[0] = 0.0
        udd_i1[0] = 0.0
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    
    i = 0
    
    while i<10 and np.max(np.abs(udd_i1-udd_i))>eps_HHT :
        
<<<<<<< HEAD
        Adh_last = Adh.copy()
        
        c_i = np.zeros((2,2))
        if np.abs(ud_i[1]) > ud_reg :
            c_i[1][1] = (Cw*m*dl_l[1]*(((L_norm/T_norm)*np.abs(ud_i[1]))**(m-1)))/M_norm
        else :
            c_i[1][1] = (Cl*dl_l[1])/M_norm
=======
        c_i = np.zeros((2,2))
        if np.abs(ud_i[1]) > ud_reg :
            c_i[1][1] = cw_norm*m*dl_norm[1]*(np.abs(ud_i[1])**(m-1))
        else :
            c_i[1][1] = cl_norm*dl_norm[1]
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
        
        ms_i = m_norm + mc*c_i + mk*k_norm
        
        ## calcul de glacier.Ff_coul2(un_i,udn_i,udn,Fc[n+1],Adh_last,n) avec glacier.F_result_typelaw2(...)
        fe_last[1] = -k_norm[1][1]*(u_i1[1] - u_last[n+1])
        fe = -np.dot(k_norm,np.transpose(u_i1))
<<<<<<< HEAD
        fp_first[0] = fp[n+1]
=======
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
        
        f_trial_i = fe + fe_last + poids_x + fp_first
        
        ## définition des états de glissement et adhérence au premier bloc de Coulomb
        if Adh[0] == False :
<<<<<<< HEAD
            if ud_i1[0]==0 or udn[0]*ud_i1[0]<0 or np.abs(f_trial_i[0]) <= p_gli: # passage de glissement à adhérence
                print('Test adhérence completed au pas :' +str(n) + ' itération '+ str(i+1))
                Adh[0] = True
        elif Adh[0] == True :
            if np.abs(f_trial_i[0]) > p_gli : # passage d'adhérence à glissement
                print('Test glissement completed au pas :' +str(n) + ' itération '+ str(i+1))
                Adh[0] = False
                if ud_i1[0] != 0:
=======
            if ud_i1[0]==0.0 or udn[0]*ud_i1[0]<0 or np.abs(f_trial_i[0]) <= p_gli: # passage de glissement à adhérence
                Adh[0] = True
        elif Adh[0] == True :
            if np.abs(f_trial_i[0]) > p_gli : # passage d'adhérence à glissement
                Adh[0] = False
                if ud_i1[0] != 0.0:
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
                    sign_last[0] = np.sign(ud_i1[0])
                else :
                    sign_last[0] = np.sign(f_trial_i[0])
        
        ## frottement de Weertman pour le deuxième bloc
<<<<<<< HEAD
        ff_weert_list = [0]
        if np.abs(ud_i1[1])>ud_reg :
            ff_weert_list.append(-(np.sign(ud_i1[1])*Cw*dl_l[1]*(((L_norm/T_norm)*np.abs(ud_i1[1]))**(m)))/F_norm)
        else :
            ff_weert_list.append(-(Cl*dl_l[1]*(L_norm/T_norm)*ud_i1[1])/F_norm)
        ff_weert = np.array(ff_weert_list)
        
        ## calcul des forces Coulomb et calcul des résultantes
        ff_alpha_list = [] # initialisation avec le frottement aux blocs de Weertman seulement
        p_alpha_t = np.zeros(2)
        
        if Adh[0] == True : 
            ff_alpha_list.append(f_trial_i[0]*sign_last[0])
            p_alpha_t[0] = 0
        if Adh[0] == False : 
            ff_alpha_list.append(-p_gli*sign_last[0])
            p_alpha_t[0] = ff_alpha[0] + f_trial_i[0]
        
        ff_alpha_list.append(ff_weert[1])
        ff_alpha = np.array(ff_alpha_list)
        
=======
        ff_weert = np.array([0,0])
        if np.abs(ud_i1[1])>ud_reg :
            ff_weert[1]=-np.sign(ud_i1[1])*cw_norm*dl_norm[1]*(np.abs(ud_i1[1])**(m))
        else :
            ff_weert[1]=-cl_norm*dl_norm[1]*ud_i1[1]
        
        ## calcul des forces Coulomb et calcul des résultantes
        ff_alpha = ff_weert # initialisation avec le frottement aux blocs de Weertman seulement
        p_alpha_t = np.zeros(2)
        
        if Adh[0] == True : 
            ff_alpha[0] = -f_trial_i[0]*sign_last[0]
            p_alpha_t[0] = 0
        if Adh[0] == False : 
            ff_alpha[0] = -p_gli*sign_last[0]
            p_alpha_t[0] = ff_alpha[0] + f_trial_i[0]
        
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
        p_alpha_t[1] = f_trial_i[1] + ff_alpha[1]
        
        u_i = u_i1
        ud_i = ud_i1
        udd_i = udd_i1
        
        # calcul de la nouvelle position avec les deux forces équilibrées mais que ce passe t'il quand décollement au pas de temps suivant P_naplha, itération, et décollement au pas de temps d'après pour Pn Pn ne change pas dans la boucle
        r_i = p_alpha*(p_n) + p_alpha1*(p_alpha_t) - np.dot(m_norm,udd_i)
        delta_uddn_i = np.dot(np.linalg.inv(ms_i),r_i)
    
        u_i1 = u_i + uiai*delta_uddn_i
        ud_i1 = ud_i + udiai*delta_uddn_i
        udd_i1 = udd_i + delta_uddn_i
        
        if Adh[0] == True :
            u_i1[0] = u_i[0]
<<<<<<< HEAD
            ud_i1[0] = 0
            udd_i1[0] = 0
=======
            ud_i1[0] = 0.0
            udd_i1[0] = 0.0
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    
        i +=1
    
    if i==10 :
        print('Alerte non convergence au pas ' + str(n))
    
<<<<<<< HEAD
    ff = ff_alpha.copy()
    p_n = p_alpha_t.copy()
    
    Adh_last = Adh.copy()
    del Adh
=======
    ff = ff_alpha
    p_n = p_alpha_t
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
    
    en = np.abs(ud_i1 - ud_i) #residual error
    
    # state vector at n+1
    u = u_i1
    ud = ud_i1
    udd = udd_i1
    
<<<<<<< HEAD
    if n%10==0 :
        u_plot.append( u )
        ud_plot.append( ud )
        udd_plot.append( udd )
        
        U_plot.append( u*L_norm )
        Ud_plot.append( ud*(L_norm/T_norm) )
        Udd_plot.append( udd*(L_norm/(T_norm*T_norm)) )
        
        ff_plot.append( ff )
        f_plot.append( p_alpha_t )
        f_trial_plot.append( f_trial_i )
        
        Ff_plot.append( ff*F_norm )
        F_plot.append( p_alpha_t*F_norm )
        F_trial_plot.append( f_trial_i*F_norm )
        Adh_list.append( Adh_last )
        Sign_list.append( sign_last ) 
        
        Niter_implicite.append(i)
        ErrResi.append(en)
=======
    u_plot.append( u )
    ud_plot.append( ud )
    udd_plot.append( udd )
    
    U_plot.append( u*L_norm )
    Ud_plot.append( ud*(L_norm/T_norm) )
    Udd_plot.append( udd*(L_norm/(T_norm*T_norm)) )
    
    ff_plot.append( ff )
    f_plot.append( p_alpha_t )
    f_trial_plot.append( f_trial_i )
    
    Ff_plot.append( ff*F_norm )
    F_plot.append( p_alpha_t*F_norm )
    F_trial_plot.append( f_trial_i*F_norm )
    Adh_list.append( Adh )
    
    Niter_implicite.append(i)
    ErrResi.append(en)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc

save_figure_path = 'C:\\users\\Paul Beguin\\Desktop\\Stage IPGP\\Modélisation Glacier Python\\Figures Python\\2020_10_19_CoulWeert_2blocs_test_CcTan'
os.chdir(save_figure_path)

import matplotlib
matplotlib.rcParams['figure.figsize'] = 35,18
matplotlib.rcParams.update({'font.size': 20})
title_font = {'size':'30'}
label_font = {'size':'30'}

<<<<<<< HEAD
T_plot = np.linspace(0,T_end,len(U_plot))

Cc_label = 'Tan_init2_dt000000314'
plt.figure(1)
plt.plot(T_plot,np.array(U_plot)*1e+3,linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([u_last[0],u_last[-1]])*L_norm*1e+3,'b--',linewidth=2)
=======
Cc_label = 'Tan'
plt.figure(1)
plt.plot(T_plot,np.array(U_plot)*1e+3,linewidth=2)
plt.plot(T_plot,u_last*L_norm*1e+3,'b--',linewidth=2)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
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
<<<<<<< HEAD
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*Ud0_last*1e+3,'y--',linewidth=2)
=======
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([Ud0_last*1e+3,Ud0_last*1e+3]),'y--',linewidth=2)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
# plt.plot([T_Fp,T_Fp],[-5,5],'r--',linewidth=2)
plt.xlabel('Temps ($s$)',label_font)
plt.ylabel('Vitesse ($mm/s$)',label_font)
plt.title("Vitesse du modèle 2 blocs avec $\mu = tan( alpha ) $",title_font)
plt.legend(['Bloc Coulomb','Bloc Weertman','$V_0$'],loc='best')
plt.grid(True)
# plt.axis([0,20,-5,5])
plt.savefig('VitesseTestCc' + Cc_label + '.png')
plt.savefig('VitesseTestCc' + Cc_label + '.pdf')

<<<<<<< HEAD
fp_plot = []
for k in range(len(T_plot)):
    if 10*k < len(T_plot) :
        fp_plot.append(fp[10*k])
    else :
        fp_plot.append(fp[-1])

plt.figure(3)
plt.plot(T_plot[1:],np.array(Ff_plot)*1e-6,linewidth=2)
plt.plot(T_plot,fp_plot*F_norm*1e-6,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),-np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.ones(2)*P_gli*1e-6,'b--',linewidth=2)
=======
plt.figure(3)
plt.plot(T_plot[1:],np.array(Ff_plot)*1e-6,linewidth=2)
plt.plot(T_plot,fp*F_norm*1e-6,'y--',linewidth=2)
# plt.plot([T_Fp,T_Fp],[-40,60],'r--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([-P_gli*1e-6,-P_gli*1e-6]),'b--',linewidth=2)
plt.plot(np.array([T_plot[0],T_plot[-1]]),np.array([P_gli*1e-6,P_gli*1e-6]),'b--',linewidth=2)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc
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

<<<<<<< HEAD
# plt.figure(5) # déplacement normalisé
# plt.plot(t_norm,u_plot)
# plt.grid(True)
# 
# plt.figure(6) # vitesse normalisée
# plt.plot(t_norm,ud_plot)
# plt.grid(True)
# 
# plt.figure(7) # accélération normalisée
# plt.plot(t_norm,udd_plot)
# plt.grid(True)
# 
# plt.figure(8) # forces de frottement
# plt.plot(t_norm[1:],ff_plot)
# plt.grid(True)
# 
# plt.figure(9) # résultantes des forces
# plt.plot(t_norm[1:],f_plot)
# plt.grid(True)
=======

plt.figure(5) # déplacement normalisé
plt.plot(t_norm,u_plot)
plt.grid(True)

plt.figure(6) # vitesse normalisée
plt.plot(t_norm,ud_plot)
plt.grid(True)

plt.figure(7) # accélération normalisée
plt.plot(t_norm,udd_plot)
plt.grid(True)

plt.figure(8) # forces de frottement
plt.plot(t_norm[1:],ff_plot)
plt.grid(True)

plt.figure(9) # résultantes des forces
plt.plot(t_norm[1:],f_plot)
plt.grid(True)
>>>>>>> d10d8e15e318b4910cabefb447bc81ad64e8cefc

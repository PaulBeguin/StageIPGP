# -*- coding: utf-8 -*-
### Modélisation SetUp1 test 2 implentation 
#reprise du code main_glacier_test8
# pour afficher le graphe 2D de la déformation u(x,t)
# avec schéma de Verlet initial
# avec implémentation de la nouvelle géométrie - glacier qui repose sur le
# avec loi de Weertman
# close all

## import des modules utiles
import numpy as np
import os
import time


## construction class Glacier pour initialisation des valeurs et fonctions
class Glacier:
    
# Paramètre de construction physique:
#   masse volumique glace : rho_glace [kg.m^-3]
#   masse volumique eau : rho_eau [kg.m^-3]
#   épaisseur du glacier : H [m]
#   pente du glacier : alpha_ground [deg]
#   accélaration de la pesanteur : g [m.s^-2]
#   coeff de friction Weertman : Cw [Pa]
#   coeff de friction Budd : Cb [Pa]
#   coeff de friction Schoof : Cs [Pa]
#   exposant de vitesse de glissement : m
#   coeff de glissement max Schoof : Cmax [Pa]
#   loi utilisée : type_law =0 Weertman, =1 Tsai, =2 sans frottement

# Lecture des fichiers utiles récupération des données
#   nom du fichier avec liste temps T : T_filename
#   nom du fichier avec liste force contact Fc : Fc_filename
#   nom du fichier avec Nt_retour indice du retournement : Nt_retour_filename
#   nom du fichier avec T_retour durée du retournement : T_retour_filename

# Construction de la géométrie du glacier et discrétisation
#   nombre de bloc : Np
#   longueur du glacier : Ltot
#   profondeur d'immergement au front du glacier : h_imm 

# Initialisation des paramètres d'optimisation
#   pas d'indice de bloc gardé pour affichage : Np_r
#   pas d'indice de temps gardé pour affichage : Nt_r
    def __init__(self,rho_glace,rho_eau,H,alpha_ground,alpha_surface,g,E,Np,Ltot,h_im):
        
        ## Paramètre de construction physique:
        # paramètres géométrique et physique élémentaires
        self.rho_glace=rho_glace #masse volumique glace
        self.rho_eau=rho_eau #masse volumique eau
        self.alpha_ground=alpha_ground #pente du glacier
        self.alpha_surface=alpha_surface
        self.g=g #accélération pesanteur
        self.E=E #module d'young de la glace

        
        ## Construction élémentaire de la géométrie du glacier et discrétisation
        # discrétisation du glacier
        self.H=H #épaisseur du glacier
        self.Ltot=Ltot #longueur du glacier
        self.Np=Np #nombre de blocs
        
        ## Setting of the glacier geometrie
        self.l_def()
        self.H_def()
        
        ## Mass matrix and mass of the blocs array
        self.M_bloc_def()
        self.M_def()
        
        ## Construction of matrxi stiffness
        self.K_def() #enregistrement utile pour calcul de Fe
        
    
    ## stiffness matrix
    def K_def(self):
        
        K_ = np.zeros((self.Np+1,self.Np+1))
        for i in range(self.Np):
            Ke = (self.E*self.H_l[i+1]/self.dl) * np.array([[1,-1],[-1,1]])
            K_[i:i+2,i:i+2] = K_[i:i+2,i:i+2] + Ke
        
        self.K = K_
        
    
    
    ## Length array for each bloc
    def l_def(self):
        
        l_init = np.linspace(0,self.Ltot,self.Np+1) #definition of the borders of the bloc
        
        dl_ = l_init[1] - l_init[0] #bloc length
        dl_l_ = np.zeros(self.Np+1)
        l_= np.zeros(self.Np+1)
        dl_l_[0] = dl_/2
        l_[0] = dl_l_[0]
        for k in range(1,self.Np):
            dl_l_[k] = dl_
            l_[k] = l_[k-1] + dl_l_[k]
        dl_l_[self.Np] = dl_/2
        l_[self.Np] = l_[self.Np-1] + dl_l_[self.Np]
        
        self.dl = dl_
        self.l = l_
        self.dl_l = dl_l_
    
    
    ## Heigth of borders 
    def H_def(self):
        
        Y_down_ = np.zeros(self.Np+2)
        Y_up_ = np.zeros(self.Np+2)
        
        Y_down_[0] = 0
        Y_up_[0] = self.H
        for k in range(1,self.Np+2):
            Y_down_[k] = np.sin(self.alpha_ground)*self.l[k-1]
            Y_up_[k] = np.sin(self.alpha_surface)*self.l[k-1] + self.H
        
        H_ = Y_up_ - Y_down_
        self.H_l = H_
        
    ## Mass for each bloc function
    def M_bloc_def(self):
        
        M_bloc_ = np.zeros(self.Np+1)
        for k in range(self.Np+1):
            M_bloc_[k] = (self.H_l[k]+self.H_l[k+1])/2 * self.rho_glace * self.dl_l[k]
        
        self.M_bloc = M_bloc_
    
    
    ## Mass matrix construction
    def M_def(self):
        
        M_ = np.diag(self.M_bloc)
        
        self.M = M_

    
    ## Matrice de réduction des déplacements de Np+1 à Np_plot+1
    def R_plot(self):
        
        R_plot_=np.zeros((self.Np_plot+1 ,self.Np+1))
        R_plot_[0][0] = 1
        R_plot_[self.Np_plot][self.Np] = 1
        for i in range(self.Np_plot-1):
            lp = self.lp_plot[i]+1
            R_plot_[i+1 ][lp]=1
        
        self.R_plot = R_plot_
    
    ## reduction matrix for sliding force
    def Rf_plot(self):
        
        Rf_plot_=np.zeros((self.Np_plot+1,self.Np+1))
        dl_mean = self.Ltot/(self.Np_plot + 1)
        Rf_plot_[0][0] = dl_mean*(2/self.dl)
        Rf_plot_[self.Np_plot][self.Np] = dl_mean*(2/self.dl)
        for i in range(self.Np_plot-1):
            lp = self.lp_plot[i]+1
            Rf_plot_[i+1][lp]= dl_mean/self.dl_l[lp]
        
        self.Rf_plot = Rf_plot_
    
    ## force elastique
    def Fe(self,U_i_):
        
        return(-np.dot(self.K,U_i_))
    
    ## force de perturbation
    def Fpertur(self,Fc_i_):
        
        # retourne la force de perturbation - avec le contact avec l'iceberg
        # i itération du temps, Fc_i_ force de contact à l'instant ti, self.Np longueur du vecteur
        
        Fpertur = np.zeros(self.Np+1)
        Fpertur[0] = Fc_i_
        
        return(Fpertur)


## Lecture de la liste T temps
def T_lecture(T_filename,path):
    
    os.chdir(path)
    T_file=open(T_filename,'r') 
    T_txt=T_file.read()
    T_file.close()
    T_txt_list=T_txt.splitlines()
    Nt=len(T_txt_list) #longueur nombre d'itération en temps
    T=[]
    for k in range(Nt):
        T.append(float(T_txt_list[k]))
        
    dt=T[1]-T[0] #pas de temps
    Ttot=T[-1] #temps de simulation
    
    return(T, dt, Ttot, Nt)
    
## Lecture de la Force de contact a la longueur Nt avec Ttot dt paramètre de T
def Fc_lecture(Fc_filename,path,dt,Ttot):
    
    # lecture de la liste 
    Ttot_fc = 300
    os.chdir(path)
    Fc_file=open(Fc_filename,'r')
    Fc_txt=Fc_file.read()
    Fc_file.close()
    Fc_txt_list=Fc_txt.splitlines()
    Fc_init=[]
    Len_Fc_load = int(Ttot*len(Fc_txt_list)/Ttot_fc)
    for k in range(Len_Fc_load):
        Fc_init.append(float(Fc_txt_list[k]))
    
    Nt_init = len(Fc_init)
    T_init = np.linspace(0,Ttot,Nt_init)
    dt_init = T_init[1]-T_init[0]
    
    # construction de la liste Fc lonugueur Nt
    T = np.arange(0,Ttot,dt)
    Nt = len(T)
    Fc = np.zeros(Nt)
    i = 0
    for k in range(Nt):
        tk = T[k]
        
        while T_init[i+1]<tk and i<Nt_init-1:
            i+=1
        
        ti = T_init[i]
        ti1 = T_init[i+1]
        Fc[k] = (ti1-tk)/dt_init *Fc_init[i] + (tk-ti)/dt_init *Fc_init[i+1]
    
    return(Fc)

## Lecture de valeur dans fichier - Nt_retour et T_retour
def Value_lecture(Value_filename,type_value,path):
    
    os.chdir(path)
    #type_value pour différencier Nt_retour de T_retour - int=0 et float=1
    Value_file = open(Value_filename,'r')
    Value_txt = Value_file.read()
    
    if type_value==0: #Nt_retour
        return(int(Value_txt.lstrip()))
    if type_value==1 : #T_retour
        return(float(Value_txt.lstrip()))
    else :
        T_max_read = Value_txt.splitlines()
        return(float(T_max_read[0]),float(T_max_read[1]))

## Construction des paramètres de réduction du temps pour affichage
def R_temps(Nt_,Nt_r_,Ttot):

    li_plot_ = np.arange(0,Nt_+Nt_r_,Nt_r_)
    Nt_plot_ = len(li_plot_)
    T_plot_ = np.linspace(0,Ttot,Nt_plot_)
    
    return(Nt_plot_,li_plot_,T_plot_)

## research the excited file names and writing of a new one
def create_file_result(glacier,results_path,theme_name):
    
    results_filename = "Results" + theme_name + "dl"+ str(int(glacier.dl))+ "Ltot" + str(int(glacier.Ltot/1000)) + "kmHHTsansFrott.npz"
    
    os.chdir(results_path)
    file_list = os.listdir(results_path)
    
    i = 2
    while results_filename in file_list :
        results_filename = "Results" + theme_name + "dl"+ str(int(glacier.dl))+ "Ltot" + str(int(glacier.Ltot/1000)) + "kmHHTsansFrottV" + str(i) +".npz"
        i +=1
    
    return(results_filename)


def create_folder_figure(glacier,save_figure_path,theme_name):
    
    os.chdir(save_figure_path)
    folder_list = os.listdir(save_figure_path)
    
    date = time.asctime(time.localtime())[0:10].replace(" ", "_")
    save_folder_name = date + "_" + theme_name
    
    i=2 
    while save_folder_name in folder_list :
        save_folder_name = date + "_" + theme_name + "_V" + str(i)
        i +=1
    
    os.makedirs(save_folder_name)
    return(save_folder_name)
    
## Fonction pour intégration temporelle et simulation
def simu_alpha(glacier,Fc,dt,Nt,alpha,eps,Nt_r):

    # définition des vecteurs de calcul et de l'état initial du système
    Un = np.zeros(glacier.Np+1) # displacement
    Udn = np.zeros(glacier.Np+1) # vitesse 
    Uddn = np.zeros(glacier.Np+1) # acceleration
    
    # arrays to return and at plotting 
    Ut_plot = [[0 for k in range(glacier.Np+1)]] # displacement
    Udt_plot = [[0 for k in range(glacier.Np+1)]] # vitesse affichée
    Uddt_plot = [[0 for k in range(glacier.Np+1)]] # acceleration
    Ftf_plot = [0] # force
    Ftsismique_plot = [0] # force de frottement comparée au frottement statique
    Ftsismique_map = [[0 for k in range(glacier.Np+1)]]
    
    Niter_implicite = []
    
    # matrix and coefficient
    mc = (1+alpha)*(0.5-alpha)*dt
    mk = (1+alpha)*0.25*((1-alpha)*dt)**2
    uiai = 0.25*((1-alpha)*dt)**2
    udiai = (0.5-alpha)*dt
    p_alpha = -alpha
    p_alpha1 = (1+alpha)
    
    M = glacier.M 
    K = glacier.K 
    
    for n in range(Nt-1):
        
        un_i = Un + dt*Udn + 0.5*(dt**2)*Uddn
        udn_i = Udn + dt*Uddn
        uddn_i = Uddn
        
        Ms_i = M + mk*K
        
        P_n = glacier.Fpertur(Fc[n]) + glacier.Fe(Un)
        P_nalpha = glacier.Fpertur(Fc[n+1]) + glacier.Fe(un_i)
        
        r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
        delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
        
        un_i1 = un_i + uiai*delta_uddn_i
        udn_i1 = udn_i + udiai*delta_uddn_i
        uddn_i1 = uddn_i + delta_uddn_i
        
        un_i1[-1] = 0
        
        i = 0
        
        while i<100 and np.max(np.abs(udn_i1 - udn_i))>eps :
            
            un_i = un_i1
            udn_i = udn_i1
            uddn_i = uddn_i1
            
            Ms_i = M + mk*K
            
            P_nalpha = glacier.Fpertur(Fc[n+1]) + glacier.Fe(un_i)
            
            r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
            
            un_i1[-1] = 0
            
            i +=1

        # state vector at n+1
        Un = un_i1
        Udn = udn_i1
        Uddn = uddn_i1
        
        Niter_implicite.append(i)
        
        # enregistrement des données pour affichage
        if n%Nt_r==0:
            print(n) 
        
        Ut_plot.append( Un )
        Udt_plot.append( Udn )
        Uddt_plot.append( Uddn )
    
        
    return(np.array(Ut_plot), np.array(Udt_plot), np.array(Uddt_plot),np.array(Niter_implicite))
    

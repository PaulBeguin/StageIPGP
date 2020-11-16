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
    def __init__(self,rho_glace,rho_eau,H_filename,alpha_ground,alpha_surface,g,E,Cw,m,Cc,type_law,ef,Np_init,Ltot,H,h_im,Xc):
        
        ## Paramètre de construction physique:
        # paramètres géométrique et physique élémentaires
        self.rho_glace=rho_glace #masse volumique glace
        self.rho_eau=rho_eau #masse volumique eau
        self.alpha_ground=alpha_ground #pente du glacier
        self.alpha_surface=alpha_surface
        self.g=g #accélération pesanteur
        self.E=E #module d'young de la glace
        
        #definition de la loi de frottement
        self.Cw=Cw
        self.Cb=Cw
        self.m=m
        self.Cc=Cc
        self.type_law=type_law
        self.ef=ef
        self.Cl=self.Cw*((self.ef)**(self.m-1))
        
        self.H = H
        self.h_im = h_im
        self.Xc = Xc
        
        ## Construction élémentaire de la géométrie du glacier et discrétisation
        # discrétisation du glacier
        self.H_filename=H_filename #file with glacier geometrie
        self.Ltot=Ltot #longueur du glacier
        self.Np_init = Np_init #nombre de blocs
        
        ## Setting of the glacier geometrie
        self.read_H_file_get_coor()
      
        ## Mass matrix and mass of the blocs array
        self.mass_def()
        self.M_bloc_def()
        self.M_def()
        

        ## Construction of matrix stiffness
        self.liaison_def()
        self.K_def() #enregistrement utile pour calcul de Fe
        
        ## calling for all the function 
        self.P_hydro_def()
        self.P_xpoids_def()
        self.P_zpoids_def()
        
        self.P_eff_def()
        self.lim_adh_def()
        
        self.U_eq_def()
        # self.Ud_eq_def()
        
        # calcul of C damping linearized matrix at steady state
        # if self.type_law==0:
        #     self.C_lin_eq = np.zeros((self.Np,self.Np))
        # if self.type_law==1:
        #     self.C_lin_eq = np.zeros((self.Np,self.Np))
    
    ## read H_filename and definition of the coor, width and hight 
    def read_H_file_get_coor(self):
        
        X_down_, Y_down_, X_up_, Y_up_ = [], [], [], []
        X_down_init, Y_down_init, X_up_init, Y_up_init = [], [], [], []
        
        file_H = open(self.H_filename,'r')
        file_H_txt = file_H.read()
        file_H.close()
        file_H_lines = file_H_txt.splitlines()
        for k in range(len(file_H_lines)):
            file_H_linesk = file_H_lines[k].split()
            X_down_init.append(float(file_H_linesk[0]))
            Y_down_init.append(float(file_H_linesk[1]))
            X_up_init.append(float(file_H_linesk[2]))
            Y_up_init.append(float(file_H_linesk[3]))
        
        dl_init = (X_down_init[0] - X_down_init[-1])/self.Np_init
        self.dl_init = dl_init
        xt = X_down_init[0]
        
        if self.type_law==1 or self.type_law==2:
            test_coul_bord = True 
            test_weert_bord = True
        if self.type_law==0 or self.type_law==3:
            test_coul_bord = False
            test_weert_bord = False
        
        if xt >= self.Ltot :
            h_im = -Y_down_init[0]
            X_down_.append(xt)
            Y_down_.append(-h_im)
            X_up_.append(xt)
            Y_up_.append(self.H-h_im)
            for k in range(self.Np_init) :
                X_down_.append(xt - (dl_init/2 + k*dl_init))
                X_up_.append(xt - (dl_init/2 + k*dl_init))
                Y_down_.append(-h_im)
                Y_up_.append(self.H-h_im)
            X_down_.append(X_down_init[-1])
            X_up_.append(X_up_init[-1])
            Y_down_.append(Y_down_init[-1])
            Y_up_.append(Y_up_init[-1])
        
        if xt <= 0:
            X_down_.append(xt)
            Y_down_.append(Y_down_init[0])
            X_up_.append(xt)
            Y_up_.append(Y_up_init[0])
            
            k_glacier = 0
            if test_coul_bord :
                self.lim_ind_coul = k_glacier
                test_coul_bord = False
            
            while (xt - (dl_init/2 + k_glacier*dl_init)) > X_down_init[-1]:
                X_down_.append(xt - (dl_init/2 + k_glacier*dl_init))
                X_up_.append(xt - (dl_init/2 + k_glacier*dl_init))
                Y_down_.append(Y_down_init[0] + np.tan(self.alpha_ground)*(dl_init/2 + k_glacier*dl_init - xt))
                Y_up_.append(Y_up_init[0] + np.tan(self.alpha_surface)*(dl_init/2 + k_glacier*dl_init - xt))
                k_glacier +=1
                
                if test_weert_bord :
                    if self.Xc > (xt - (dl_init/2 + k_glacier*dl_init)) :
                        X_down_.append(self.Xc)
                        X_up_.append(self.Xc)
                        Y_down_.append(Y_down_init[0] + np.tan(self.alpha_ground)*(-self.Xc))
                        Y_up_.append(Y_up_init[0] + np.tan(self.alpha_surface)*(-self.Xc))
                        
                        test_weert_bord = False
                    self.lim_ind_weert = k_glacier
            
            X_down_.append(X_down_init[-1])
            X_up_.append(X_up_init[-1])
            Y_down_.append(Y_down_init[-1])
            Y_up_.append(Y_up_init[-1])
        
        if xt<self.Ltot and xt>0 :
            X_down_.append(xt)
            X_up_.append(xt)
            Y_down_.append(Y_down_init[0])
            Y_up_.append(Y_up_init[0])
            
            k_ton = 0
            while (xt - (dl_init/2 + k_ton*dl_init)) > 0:
                X_down_.append(xt - (dl_init/2 + k_ton*dl_init))
                X_up_.append(xt - (dl_init/2 + k_ton*dl_init))
                Y_down_.append(Y_down_init[0])
                Y_up_.append(Y_up_init[0])
                k_ton +=1
            
            k_glacier = k_ton
            if test_coul_bord :
                self.lim_ind_coul = k_glacier
                test_coul_bord = False
            
            if (xt - (dl_init/2 + k_glacier*dl_init)) < 0 :
                X_down_.append(0)
                X_up_.append(0)
                Y_down_.append(Y_down_init[1])
                Y_up_.append(Y_up_init[1])
            
            while (xt - (dl_init/2 + k_glacier*dl_init)) > X_down_init[-1]:
                X_down_.append(xt - (dl_init/2 + k_glacier*dl_init))
                X_up_.append(xt - (dl_init/2 + k_glacier*dl_init))
                Y_down_.append(Y_down_init[1] + np.tan(self.alpha_ground)*(dl_init/2 + k_glacier*dl_init - xt))
                Y_up_.append(Y_up_init[1] + np.tan(self.alpha_surface)*(dl_init/2 + k_glacier*dl_init - xt))
                k_glacier +=1
                
                if test_weert_bord :
                    if self.Xc > (xt - (dl_init/2 + k_glacier*dl_init)) :
                        X_down_.append(self.Xc)
                        X_up_.append(self.Xc)
                        Y_down_.append(Y_down_init[1] + np.tan(self.alpha_ground)*(-self.Xc))
                        Y_up_.append(Y_up_init[1] + np.tan(self.alpha_surface)*(-self.Xc))
                        
                        test_weert_bord = False
                    self.lim_ind_weert = k_glacier
            
            X_down_.append(X_down_init[-1])
            X_up_.append(X_up_init[-1])
            Y_down_.append(Y_down_init[-1])
            Y_up_.append(Y_up_init[-1])
        
        self.X_down = np.array(X_down_)
        self.Y_down = np.array(Y_down_)
        self.X_up = np.array(X_up_)
        self.Y_up = np.array(Y_up_)
        
        self.Np = len(X_down_) - 1
        
        H_l_ = np.zeros(len(self.Y_up))
        for k in range(len(self.Y_up)):
            H_l_[k] = self.Y_up[k] - self.Y_down[k]
        
        self.H_l = H_l_
    
    def liaison_def(self):
        
        l_liaison_ = np.zeros(self.Np-1)
        H_liaison_ = np.zeros(self.Np-1)
        
        l_liaison_[0] = np.sqrt((self.X_down[0]-(self.X_down[1]+self.X_down[2])/2)**2 + (self.Y_down[0]-(self.Y_down[1]+self.Y_down[2])/2)**2 )
        H_liaison_[0] = self.Y_up[1] - self.Y_down[1]
        
        if (self.Np-3) > 0:
            for k in range(1,self.Np-2):
                
                l_liaison_[k] = np.sqrt(((self.X_down[k]-self.X_down[k+2])/2)**2 + ((self.Y_down[k]-self.Y_down[k+2])/2)**2 )
                H_liaison_[k] = self.Y_up[k+2] - self.Y_down[k+2]
        
        l_liaison_[-1] = (np.sqrt((self.X_down[-1]-(self.X_down[-2]+self.X_down[-3])/2)**2 + (self.Y_down[-1]-(self.Y_down[-2]+self.Y_down[-3])/2)**2 ))
        H_liaison_[-1] = (self.Y_up[-2] - self.Y_down[-2])
        
        self.l_liaison = l_liaison_
        self.H_liaison = H_liaison_
    
    
    def mass_def(self):
        
        l_mass_ = np.zeros(self.Np)
        H_mass_ = np.zeros(self.Np)
        for k in range(self.Np):
            l_mass_[k] = (self.X_down[k] - self.X_down[k+1])
            H_mass_[k] = ((self.H_l[k]+self.H_l[k+1])/2)
        
        self.l_mass = l_mass_
        self.H_mass = H_mass_
        self.dl_l = l_mass_ 
    
    
    ## stiffness matrix
    def K_def(self):
        
        K_ = np.zeros((self.Np, self.Np))
        
        for i in range(self.Np-1):
            Ke = (self.E*self.H_liaison[i]/self.l_liaison[i]) * np.array([[1,-1],[-1,1]])
            K_[i:i+2,i:i+2] = K_[i:i+2,i:i+2] + Ke
        
        self.K = K_
    
    
    ## Mass for each bloc function
    def M_bloc_def(self):
        
        M_bloc_ = np.zeros(self.Np)
        
        for k in range(self.Np):
            M_bloc_[k] = (self.H_l[k]+self.H_l[k+1])/2 * self.rho_glace * self.dl_l[k]
        
        self.M_bloc = M_bloc_
    
    
    ## Mass matrix construction
    def M_def(self):
        
        M_ = np.diag(self.M_bloc)
        
        self.M = M_
    
    
    ## pression hydrostatique sur les blocs
    def P_hydro_def(self):
        
        P_hydro_ = np.zeros(self.Np)
        for k in range(self.Np):
            h0 = -self.Y_down[k]
            h1 = -self.Y_down[k+1]
            if h0>0 :
                if h1>=0 :
                    P_hydro_[k] = self.rho_eau * self.g * self.dl_l[k] * (h0+h1)/2
                elif h1<0 :
                    dlk = self.dl_l[k] * h1/(h1-h0)
                    P_hydro_[k] = self.rho_eau * self.g * h1 * dlk/2
        
        self.P_hydro = P_hydro_
    
    ## pression poids composante normale au sol
    def P_zpoids_def(self):
        
        P_zpoids_ = np.zeros(self.Np)
        for k in range(self.Np):
            if self.X_down[k] <= 0 : 
                P_zpoids_[k] = -np.cos(self.alpha_ground) * self.g * self.M_bloc[k]
            else : 
                P_zpoids_[k] = -self.g * self.M_bloc[k]
        
        self.P_zpoids = P_zpoids_
    
    ## pression poids composante tangentielle au sol - et valeur de frottement sur tout le glacier
    def P_xpoids_def(self):
        
        P_xpoids_ = np.zeros(self.Np)
        for k in range(self.Np):
            if self.X_down[k] <= 0 : 
                P_xpoids_[k] = np.sin(self.alpha_ground) * self.g * self.M_bloc[k]

        self.F_stat = np.sum(P_xpoids_ )
        
        self.P_xpoids = P_xpoids_ 

    ## effective pressure 
    def P_eff_def(self):

        P_eff_ = self.negat(self.P_hydro + self.P_zpoids) 
        
        N_eff_ = np.zeros(self.Np)
        for k in range(self.Np):
            if self.P_zpoids[k] != 0 :
                N_eff_[k] = np.abs(P_eff_[k]/self.P_zpoids[k])
        
        self.P_eff = P_eff_
        self.N_eff = N_eff_
    
    ## positive part
    def posit(self,A):
        
        A_po = (np.abs(A) + A)/2
        
        return(A_po)
    
    
    ## negative part
    def negat(self,A):
        
        A_neg = (np.abs(A) - A)/2
        
        return(A_neg)
    
    
    ## strains lim of adherence on each block
    def lim_adh_def(self):
        
        self.F_coulomb_lim = self.Cc*self.P_eff
    
    ## initialisation of steady deformation state
    def U_eq_def(self):
        
        if self.type_law==0 :
            u_eq_ = np.zeros(self.Np)
            self.u_eq = u_eq_
        
        if self.type_law==1 :
            # u_eq_ = np.zeros(self.Np)
            ind1 = 0
            ind2 = 0
            for k in range(self.Np):
                if self.X_down[k]> 0:
                    ind1 = k
                if self.X_down[k]>=self.Xc :
                    ind2 = k
            # K_coul = self.K[ind1:ind2,ind1:ind2]
            F_eq_ = np.zeros(self.Np)
            F_eq_coul = (self.P_xpoids[ind1:ind2] - self.F_coulomb_lim[ind1:ind2])
            F_eq_[ind1:ind2] = F_eq_coul
            # U_eq_coul = np.dot(np.linalg.inv(K_coul),F_eq_coul)
            # 
            # if ind1 >0:
            #     u_eq_[:ind1] = U_eq_coul[0]
            # u_eq_[ind1:ind2] = U_eq_coul
            # u_eq_[ind2:] = U_eq_coul[-1]
            
            u_eq_ = np.zeros(self.Np)
            u_eq_[:-1] = np.dot(np.linalg.inv(self.K[:-1,:-1]),F_eq_[:-1])
            self.U_eq = u_eq_
    
    ## initialisation of steady sliding state
    def Ud_eq_def(self):
        
        P_eff_ = self.P_eff
        P_xpoids_ = self.P_xpoids
        
        Ud_eq_ = np.zeros(self.Np)
        if self.type_law==0:
            ind_eq_ton = self.Np-1
            test_eq_ton = True
            for k in range(self.Np):
                if self.X_down[k]<=0:
                    if test_eq_ton :
                        ind_eq_ton = k
                        test_eq_ton = False
                    Ud_eq_[k] = (P_xpoids_[k]/(self.Cw*self.dl_l[k]))**(1/self.m)
            for k in range(ind_eq_ton):
                Ud_eq_[k] = Ud_eq_[ind_eq_ton]
            
        if self.type_law==1:
            ind_eq_ton = self.Np-1
            test_eq_ton = True
            for k in range(self.Np):
                if self.X_down[k]<= self.Xc:
                    if test_eq_ton :
                        ind_eq_ton = k
                        test_eq_ton = False
                    Ud_eq_[k] = (P_xpoids_[k]/(self.Cw*self.dl_l[k]))**(1/self.m)
            for k in range(ind_eq_ton):
                Ud_eq_[k] = Ud_eq_[ind_eq_ton]
        
        self.Ud_eq = Ud_eq_
 
    ## Force élastique
    def Fe(self,U):
        
        # retourne la force élastique aux blocs
        # K=self.K , U=np.transpose(U), multiplication avec np.dot
        
        return(-np.dot(self.K,np.transpose(U)))
    
    ## sliding force with Power law
    def Ff_weert_nonreg(self,Ud):
        
        Ff_ = np.zeros(self.Np)
        
        for k in range(self.Np):
            Ff_[k]=-np.sign(Ud[k])*self.Cw*self.dl_l[k]*(np.abs(Ud[k])**(self.m))
            
        
        return(Ff_)
    
    ## sliding force with Power law
    def Ff_weert(self,Ud):
        
        Ff_ = np.zeros(self.Np)
        
        if self.type_law==0 or self.type_law==2:
            for k in range(self.Np):
                if self.X_down[k] <= 0 :
                    if np.abs(Ud[k])>self.ef :
                        Ff_[k]=-np.sign(Ud[k])*self.Cw*self.dl_l[k]*(np.abs(Ud[k])**(self.m))
                    else :
                        Ff_[k]=-self.Cl*self.dl_l[k]*Ud[k]
            
        if self.type_law==1:
            for k in range(self.lim_ind_weert,self.Np):
                if self.X_down[k] <= self.Xc :
                    if np.abs(Ud[k])>self.ef :
                        Ff_[k]=-np.sign(Ud[k])*self.Cw*self.dl_l[k]*(np.abs(Ud[k])**(self.m))
                    else :
                        Ff_[k]=-self.Cl*self.dl_l[k]*Ud[k]
        
        return(Ff_)
    
    
    
    ## sign function made for a strain array
    def sign(self,A):
        
        sign_result_ = np.zeros(len(A))
        for k in range(len(A)):
            if A[k] < 0 :
                sign_result_[k] = -1
            if A[k] > 0 :
                sign_result_[k] = 1
        
        return(sign_result_)
    
    
    ## sliding force with Coulomb law

    def Ff_coul(self,U_,U_adh_,DU_adh_,udi1,udi,eps_):
        
        F_trial_ = -eps_*(U_-(U_adh_+DU_adh_))
        # print(F_trial_)
        Ff_coul_ = np.zeros(self.Np)
        
        DDU_adh_r = DU_adh_.copy()
        adh_state = [False for k in range(glacier.Np)]
        
        for k in range(self.Np):
            if (self.X_down[k]<=0) and (self.X_down[k]>self.Xc):
                if (np.abs(F_trial_[k])<= self.F_coulomb_lim[k] and udi1[k]==0 ) or udi[k]*udi1[k] < 0 :
                    Ff_coul_[k] = F_trial_[k]
                else :
                    Ff_coul_[k] = np.sign(F_trial_[k])*np.sign(udi1[k])
                    DDU_adh_r[k] = - np.sign(F_trial_[k])*(np.abs(F_trial_[k])-self.F_coulomb_lim[k])/eps_
        
        return(Ff_coul_,DDU_adh_r,adh_state)
    
    ## sliding force with Coulomb law sans pénalisation

    def Ff_coul2(self,U_,Ud_,Ud_last_,Fc_,Adh_last,sign_last):
        
        F_trial_ = self.Fe(U_) + self.P_xpoids + self.F_pertur(Fc_)
        
        Adh_r_ = Adh_last.copy()
        sign_r_ = sign_last.copy()
        
        for k in range(self.lim_ind_coul,self.lim_ind_weert+1):
            if Adh_last[k] == False :
                if (Ud_[k]==0) or Ud_[k]*Ud_last_[k]<0 :
                    Adh_r_[k] = True
            elif Adh_last[k] == True :
                if np.abs(self.F_coulomb_lim[k]) < np.abs(F_trial_[k]) :
                    Adh_r_[k] = False
                    if Ud_[k] != 0:
                        sign_r_[k] = np.sign(Ud_[k])
                    else :
                        sign_r_[k] = np.sign(F_trial_[k])
        
        Ff_coul_ = np.zeros(self.Np)
        
        for k in range(self.lim_ind_coul,self.lim_ind_weert+1):
            if Adh_r_[k] == True :
                Ff_coul_[k] = - np.abs(F_trial_[k])*sign_r_[k]
            elif Adh_r_[k] == False :
                Ff_coul_[k] = - np.abs(self.F_coulomb_lim[k])*sign_r_[k]
        
        return(Ff_coul_,Adh_r_,sign_r_) 
    
    
    ## budd law - weertman with N_eff influence
    def Ff_budd(self,Ud_):
        
        Ff_ = np.zeros(self.Np)
        for k in range(self.Np):
            if self.X_down[k] <= 0:
                if np.abs(Ud_[k])>self.ef :
                    Ff_[k]=-np.sign(Ud_[k])*self.Cw*self.dl_l[k]*(np.abs(Ud_[k])**(self.m))*self.N_eff[k]
                else :
                    Ff_[k]=-self.Cl*self.dl_l[k]*Ud_[k]*self.N_eff[k]
        
        return(Ff_)
    
    ## strain of perturbation resulting from iceberg impact on the first bloc
    def F_pertur(self,Fc_i_):
        
        # retourne la force de perturbation - avec le contact avec l'iceberg
        # i itération du temps, Fc_i_ force de contact à l'instant ti, self.Np longueur du vecteur
        
        Fpertur = np.zeros(self.Np)
        Fpertur[0] = Fc_i_
        
        return(Fpertur)

    ## total strain for every type law
    def F_result_typelaw_nonreg(self,U_,Ud_,Fc_):
        
        F_result_ = self.P_xpoids + self.F_pertur(Fc_) + self.Fe(U_) + self.Ff_weert_nonreg(Ud_)
        # F_result_[-1] = 0
        Ff_ = self.Ff_weert(Ud_)
        
        return(F_result_,Ff_)
        
    ## total strain for every type law
    def F_result_typelaw0(self,U_,Ud_,Fc_):
        
        F_result_ = self.P_xpoids + self.F_pertur(Fc_) + self.Fe(U_) + self.Ff_weert(Ud_)
        # F_result_[-1] = 0
        Ff_ = self.Ff_weert(Ud_)
        
        return(F_result_,Ff_)
    
    
    def F_result_typelaw1(self,U_,Ud_,Ud_last_,Fc_,Adh_last,sign_last):
        
        Ff_coul_,Adh_r_,sign_r = self.Ff_coul2(U_,Ud_,Ud_last_,Fc_,Adh_last,sign_last)
        Ff_weert_ = self.Ff_weert(Ud_)
        Ff_ = Ff_weert_ + Ff_coul_
        
        Finter_ = self.P_xpoids + self.F_pertur(Fc_) + self.Fe(U_) 
        F_result_ = Finter_ + Ff_
        
        for k in range(self.lim_ind_coul,self.lim_ind_weert+1):
            if Adh_r_[k] == True :
                Ff_[k] = -Finter_[k]*sign_r[k]
                F_result_[k] = 0
            # if np.abs(Ff_[k]) > self.F_coulomb_lim[k]:
            #     Ff_[k] = -self.F_coulomb_lim[k]*np.sign(Ff_[k])
            #     F_result_[k] = Finter_[k] + Ff_[k]
            #     Adh_r_[k] = False
        
        return(F_result_,Ff_,Adh_r_,sign_r)
    
    def F_result_typelaw2(self,U_,Ud_,Ud_last_,sign_last,Fc_,Adh_last,n):
        
        Ff_coul_,sign_r_,Adh_r_ = self.Ff_coul2(U_,Ud_,Ud_last_,sign_last,Fc_,Adh_last,n)
        Ff_weert_ = self.Ff_weert(Ud_)
        Type_law_list = [0 for k in range(self.Np)]
        Ff_ = np.zeros(self.Np)
        
        Ff_inter = np.array([[Ff_weert_[k],Ff_coul_[k]] for k in range(self.lim_ind_coul,self.lim_ind_weert+1)])
        for k in range(self.lim_ind_coul,self.lim_ind_weert+1):
            Ffk = np.min(np.abs(Ff_inter),axis=k)
            if Ffk == np.abs(Ff_inter[k][1]) :
                Type_law_list[k] = 1
            else :
                Type_law_list[k] = 0
            Ff_[k] = Ffk
        
        F_result_ = self.P_xpoids + self.F_pertur(Fc_) + self.Fe(U_)  + Ff_
        
        return(F_result_,Ff_,sign_r_,Adh_r_,Type_law_list)
    
    ## total strain for every type law
    def F_result_typelaw3(self,U_,Ud_,Fc_):
        
        F_result_ = self.P_xpoids + self.F_pertur(Fc_) + self.Fe(U_) + self.Ff_budd(Ud_)
        # F_result_[-1] = 0
        Ff_ = self.Ff_budd(Ud_)
        
        return(F_result_,Ff_)
    
## Lecture de la Force de contact a la longueur Nt avec Ttot dt paramètre de T
def Fc_lecture(Fc_filename,path,T,T_start):
    
    # lecture de la liste 
    Ttot = T[-1]
    Ttot_fc = 300
    os.chdir(path)
    Fc_file=open(Fc_filename,'r')
    Fc_txt=Fc_file.read()
    Fc_file.close()
    Fc_txt_list=Fc_txt.splitlines()
    Fc_init=[]
    
    if Ttot<=Ttot_fc:
        Len_Fc_load = int(Ttot*len(Fc_txt_list)/Ttot_fc)
    else :
        Len_Fc_load = len(Fc_txt_list)
    
    for k in range(Len_Fc_load):
        Fc_init.append(float(Fc_txt_list[k]))
    
    Nt_init = len(Fc_init)
    if Ttot<=Ttot_fc:
        T_init = np.linspace(0,Ttot,Nt_init)
    else :
        T_init = np.linspace(0,Ttot_fc,Nt_init)
    dt_init = T_init[1]-T_init[0]
    
    
    # construction de la liste Fc lonugueur Nt
    Nt = len(T)
    Fc = np.zeros(Nt)
    i = 0
    for k in range(Nt):
        tk = T[k]
        if tk>=T_start and (tk-T_start)<=T_init[-1]:
            tk_offset = tk-T_start
            
            while T_init[i+1]<tk_offset and i<Nt_init-1:
                i+=1
            
            ti = T_init[i]
            ti1 = T_init[i+1]
            Fc[k] = (ti1-tk_offset)/dt_init *Fc_init[i] + (tk_offset-ti)/dt_init *Fc_init[i+1]
    
    return(Fc)


def Write_H_file_linear(alpha_ground,alpha_surface,Ltot,H0,Nh,H_filename):
    
    l_H = np.zeros(Nh+2)
    dl = Ltot/Nh
    l_H[1] = dl/2
    for k in range(1,Nh):
        l_H[k+1] = l_H[k] + dl
    l_H[-1] = l_H[Nh] + dl/2
    
    H_down = np.sin(alpha_ground)*l_H
    H_up = np.sin(alpha_surface)*l_H + H0
    
    H_txt = ''
    for k in range(len(H_down)):
        H_txt += "%.2f"%(l_H[k]) + '  ' + "%.2f" %(H_down[k]) + '  ' + "%.2f" %(H_up[k]) + '\n'
    
    file_H = open(H_filename,'w')
    file_H.write(H_txt)
    file_H.close()
    
    
def Write_coord_file_linear(alpha_ground, alpha_surface, xt, Ltot, H, h_im, H_filename):
    
    X_down = []
    Y_down = []
    X_up = []
    Y_up = []
    
    if xt>0 : #there is a tongue part
        X_down.append(xt)
        Y_down.append(-h_im)
        X_up.append(xt)
        Y_up.append(H - h_im)
        
        if Ltot <= xt : #the length isn't long enough to describe the ground support of the glacier
            X_down.append(xt-Ltot)
            Y_down.append(-h_im)
            X_up.append(xt-Ltot)
            Y_up.append(H-h_im)
            print('Caution tongue not long enough to reach the calving point')
        else :
            X_down.append(0)
            Y_down.append(-h_im)
            X_up.append(0)
            Y_up.append(H-h_im)
            
            X_down.append(-np.cos(alpha_ground)*(Ltot-xt))
            Y_down.append(-h_im + np.sin(alpha_ground)*(Ltot-xt))
            X_up.append(-np.cos(alpha_surface)*(Ltot-xt))
            Y_up.append(H-h_im + np.sin(alpha_surface)*(Ltot-xt))
    else:
        X_down.append(xt)
        Y_down.append(-h_im + np.sin(alpha_ground)*(-xt))
        X_up.append(xt)
        Y_up.append(H - h_im + np.sin(alpha_surface)*(-xt))
        
        X_down.append(xt-np.cos(alpha_ground)*Ltot)
        Y_down.append(-h_im + np.sin(alpha_ground)*(Ltot-xt))
        X_up.append(xt-np.cos(alpha_surface)*Ltot)
        Y_up.append(H-h_im + np.sin(alpha_surface)*(Ltot-xt))
        
    H_txt = ""
    for k in range(len(X_down)):
        H_txt += "%.2f"%(X_down[k]) + '  ' + "%.2f"%(Y_down[k]) + '  ' + "%.2f"%(X_up[k]) + '  ' + "%.2f"%(Y_up[k]) + '\n'
    
    file_H = open(H_filename,'w')
    file_H.write(H_txt)
    file_H.close()


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
    
    results_filename = "Results" + theme_name + "dl"+ str(int(glacier.dl_init))+ "Ltot" + str(int(glacier.Ltot/1000)) + "kmHHTtypelaw"+str(glacier.type_law)+".npz"
    
    os.chdir(results_path)
    file_list = os.listdir(results_path)
    
    i = 2
    while results_filename in file_list :
        results_filename = "Results" + theme_name + "dl"+ str(int(glacier.dl_init))+ "Ltot" + str(int(glacier.Ltot/1000)) + "kmHHTtypelaw"+str(glacier.type_law) + "V" + str(i) +".npz"
        i +=1
    
    return(results_filename)

## get the response for the perturbation only
def Get_pertur_result(results_path, results_filename_statique, results_filename_brut, results_filename):
    
    file_results_statique = np.load(results_path + '\\' + results_filename_statique)
    Ut0 = file_results_statique['Ut']
    Utd0 = file_results_statique['Utd']
    Utdd0 = file_results_statique['Utdd']
    Ftf0 = file_results_statique['Ftf']
    Ftsismique0 = file_results_statique['Ftsismique']
    Ftsismique_map0 = file_results_statique['Ftsismique_map']
    Niter_implicite0 = file_results_statique['Niter_implicite']
    file_results_statique.close()
    
    file_results_brut = np.load(results_path + '\\' + results_filename_brut)
    Ut = file_results_brut['Ut']
    Utd = file_results_brut['Utd']
    Utdd = file_results_brut['Utdd']
    Ftf = file_results_brut['Ftf']
    Ftsismique = file_results_brut['Ftsismique']
    Ftsismique_map = file_results_brut['Ftsismique_map']
    Niter_implicite = file_results_brut['Niter_implicite']
    file_results_brut.close()
    
    Utp = Ut - Ut0
    Utdp = Utd - Utd0
    Utddp = Utdd - Utdd0
    Ftfp = Ftf - Ftf0
    Ftsismiquep = Ftsismique - Ftsismique0
    Ftsismique_mapp = Ftsismique_map - Ftsismique_map0
    Niter_implicitep = Niter_implicite - Niter_implicite0
    
    np.savez(results_path + '\\' + results_filename,Ut=Utp,Utd=Utdp,Utdd=Utddp,Ftf=Ftfp,Ftsismique=Ftsismiquep,Ftsismique_map=Ftsismique_mapp,Niter_implicite=Niter_implicitep)


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
def simu(glacier,Fc,dt,Nt,Nt_r):
    
    # définition des vecteurs de calcul et de l'état initial du système
    Ui = np.zeros(glacier.Np) # déplacement - nul
    Udi = glacier.Ud_eq # vitesse - prise à l'équilibre
    Fi = glacier.Fpertur(Fc[0]) # chargement - Fc[0]
    
    # définition du système optimisé - retourné et affiché
    Ut_plot = [[0 for k in range(glacier.Np)]] # displacement
    Udt_plot = [[0 for k in range(glacier.Np)]] # vitesse affichée
    Uddt_plot = [[0 for k in range(glacier.Np)]] # acceleration
    Ftf_plot = [0] # force
    Ftsismique_plot = [0] # force de frottement comparée au frottement statique
    Ftsismique_map = [[0 for k in range(glacier.Np)]]
    
    # matrice du problème
    M = glacier.M #masse
    M_inv = np.linalg.inv(M) #calcul de l'inverse de masse
    
    # Boucle temporelle - Nt-1 itérations pour obtenir longueur finale de Ui,Fi et Udi Nt
    print(" Début d'intégration temporelle \n Rang enregistré: ") #message pour suivi du calcul
    for i in range(Nt-1):
        
        # intégration de Verlet
        Ui = Ui + dt*(Udi + (dt/2)*np.dot(M_inv,Fi))
        Fi1 = glacier.F_result_typelaw(Ui,Udi,Fc[i])
        Udi = Udi + (dt/2)*np.dot(M_inv,(Fi + Fi1))

        Fi = Fi1
        # Udi[-1] = glacier.Ud_eq[-1]
        
        # enregistrement des données pour affichage
        if i%Nt_r==0:
            print(i) 
        
        Ut_plot.append( Ui )
        Udt_plot.append( Udi )
        Uddt_plot.append( Uddi )
        Ftf_plot.append( np.sum(glacier.Ff_law0(Udi)) )
        Ftsismique_plot.append( np.sum(glacier.Ff(Udi)) + glacier.F_stat )
        Ftsismique_map.append(glacier.Ff_law0(Udi) + glacier.P_xpoids)
        
    return(np.array(Ut_plot), np.array(Udt_plot), np.array(Uddt_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map),np.array(Niter_implicite))
    
## Fonction pour intégration temporelle et simulation
def simu_alpha(glacier,Ut_init,Utd_init,Utdd_init,Fc,T,alpha,eps,eps_stiff,Nt_r):

    # information reading from the time array
    Nt = len(T)
    dt = T[1] - T[0]
    T_plot = [0]
    
    # définition des vecteurs de calcul et de l'état initial du système
    Un = Ut_init # displacement
    Udn = Utd_init # vitesse 
    Uddn = Utdd_init # acceleration
    
    # definition sticking state in displacement
    Un_adh = Ut_init # initial state
    
    # arrays to return and at plotting 
    Ut_plot = [[Un[k] for k in range(glacier.Np)]] # displacement
    Udt_plot = [[Udn[k] for k in range(glacier.Np)]] # vitesse affichée
    Uddt_plot = [[Uddn[k] for k in range(glacier.Np)]] # acceleration
    Ftf_plot = [0] # force
    Ftsismique_plot = [0] # force de frottement comparée au frottement statique
    Ftsismique_map = [[0 for k in range(glacier.Np)]]
    
    Niter_implicite = []
    ErrResi = []
    
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
        
        dt = T[n+1] - T[n] #variable step time
        
        un_i = Un + dt*Udn + 0.5*(dt**2)*Uddn
        udn_i = Udn + dt*Uddn
        uddn_i = Uddn
        
        if glacier.type_law==0:
            Cn_i = glacier.C_lin_def(udn_i)
        if glacier.type_law==1:
            Cn_i = np.zeros((glacier.Np,glacier.Np))
        
        Ms_i = M + mc*Cn_i + mk*K
        
        if glacier.type_law==0:
            P_n,Ff_n = glacier.F_result_typelaw0(Un,Udn,Fc[n]) 
            P_nalpha,Ff_nalpha = glacier.F_result_typelaw0(un_i,udn_i,Fc[n+1])
            
            r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
            
        if glacier.type_law==1:
            P_n,Ff_n = glacier.F_result_typelaw1(Un,Udn,Fc[n]) 
            P_nalpha_t,Ff_nalpha_t = glacier.F_result_typelaw1(un_i,udn_i,Fc[n+1])
            
            r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha_t) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
            
            Ff_nalpha_coul,Un_adh = glacier.F_coulomb(un_i1,Un_adh,eps_stiff)
            P_nalpha = Ff_nalpha_coul + P_nalpha_t
            
            r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
        # udn_i1[-1] = glacier.Ud_eq[-1]
        
        i = 0
        
        while i<100 and np.max(np.abs(udn_i1 - udn_i))>eps :
            
            un_i = un_i1
            udn_i = udn_i1
            uddn_i = uddn_i1
            
            if glacier.type_law==0:
                P_nalpha,Ff_nalpha = glacier.F_result_typelaw0(un_i,udn_i,Fc[n+1])
                
                r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
                delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
                
                un_i1 = un_i + uiai*delta_uddn_i
                udn_i1 = udn_i + udiai*delta_uddn_i
                uddn_i1 = uddn_i + delta_uddn_i
                
            if glacier.type_law==1:
                P_nalpha_t,Ff_nalpha_t = glacier.F_result_typelaw1(un_i,udn_i,Fc[n+1])
                
                r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha_t) - np.dot(M,uddn_i)
                delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
                
                un_i1 = un_i + uiai*delta_uddn_i
                udn_i1 = udn_i + udiai*delta_uddn_i
                uddn_i1 = uddn_i + delta_uddn_i
                
                Ff_nalpha_coul,Un_adh = glacier.F_coulomb(un_i1,Un_adh,eps_stiff)
                P_nalpha = Ff_nalpha_coul + P_nalpha_t
                
                r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
                delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
                
                un_i1 = un_i + uiai*delta_uddn_i
                udn_i1 = udn_i + udiai*delta_uddn_i
                uddn_i1 = uddn_i + delta_uddn_i
            
            # udn_i1[-1] = glacier.Ud_eq[-1]
            
            i +=1
        
        en = np.max(np.abs(udn_i1 - udn_i)) #residual error
        
        # state vector at n+1
        Un = un_i1
        Udn = udn_i1
        Uddn = uddn_i1
        
        # enregistrement des données pour affichage
        if n%Nt_r==0:
            if n%1000 == 0:
                print(n) 
            T_plot.append(T[n])
            Ut_plot.append( Un )
            Udt_plot.append( Udn )
            Uddt_plot.append( Uddn )
            
            Ftf_plot.append( np.sum(P_nalpha) )
            if glacier.type_law==0:
                Ftsismique_plot.append(  np.sum(glacier.Ff_weert(Udn)) )
                Ftsismique_map.append( glacier.Ff_weert(Udn) )
            if glacier.type_law==1:
                Ftsismique_plot.append(  np.sum(Ff_nalpha_coul) )
                Ftsismique_map.append( Ff_nalpha_coul )
            Niter_implicite.append(i)
            ErrResi.append(en)
        
    return(np.array(T_plot),np.array(Ut_plot), np.array(Udt_plot), np.array(Uddt_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map),np.array(Niter_implicite),np.array(ErrResi))


def simu_alpha_pen(glacier,Ut_init,Utd_init,Utdd_init,Fc,T,alpha_HHT,eps_HHT,eps_Pen,Nt_r):

    # T_plot0,Ut0,Utd0,Utdd0,Ftf0,Ftsismique0,Ftsismique_map0,Niter_implicite0,ErrResi0 = SeismeGlacier.simu_alpha_pen(glacier,Ut0_init,Utd0_init,Utdd0_init,Fc0,T,alpha_HHT,eps_HHT,eps_Pen,Nt_r)
    # information reading from the time array
    dt = T[1] - T[0]
    T_plot = [0]
    Nt = len(T)
    
    # définition des vecteurs de calcul et de l'état initial du système
    Un = Ut_init # displacement
    Udn = Utd_init # vitesse 
    Uddn = Utdd_init # acceleration
    
    # definition sticking state in displacement
    Un_adh = Ut_init - glacier.P_xpoids/eps_Pen # initial state
    
    # arrays to return and at plotting 
    Ut_plot = [[Un[k] for k in range(glacier.Np)]] # displacement
    Utd_plot = [[Udn[k] for k in range(glacier.Np)]] # vitesse affichée
    Utdd_plot = [[Uddn[k] for k in range(glacier.Np)]] # acceleration
    
    Ut_adh_plot = [[Un_adh[k] for k in range(glacier.Np)]]
    Ftf_plot = [0] # force
    Ftsismique_plot = [0] # force de frottement comparée au frottement statique
    Ftsismique_map = [[0 for k in range(glacier.Np)]]
    
    Niter_implicite = []
    ErrResi = []
    
    # matrix and coefficient
    mc = (1+alpha_HHT)*(0.5-alpha_HHT)*dt
    mk = (1+alpha_HHT)*0.25*((1-alpha_HHT)*dt)**2
    uiai = 0.25*((1-alpha_HHT)*dt)**2
    udiai = (0.5-alpha_HHT)*dt
    p_alpha = -alpha_HHT
    p_alpha1 = (1+alpha_HHT)
    
    M = glacier.M
    K = glacier.K
    
    for n in range(Nt-1):
        
        un_i = Un + dt*Udn + 0.5*(dt**2)*Uddn
        udn_i = Udn + dt*Uddn
        uddn_i = Uddn
    
        Cn_i = np.zeros((glacier.Np,glacier.Np))
        for k in range(glacier.lim_ind_weert+1,glacier.Np):
            if np.abs(udn_i[k])>=glacier.ef :
                Cn_i[k][k] = glacier.m*glacier.Cw*glacier.dl_l[k]*(np.abs(udn_i[k])**(glacier.m - 1))
            else :
                Cn_i[k][k] = glacier.Cl*glacier.dl_l[k]
    
        
        Ms_i = M + mc*Cn_i + mk*K
        
        P_n,Ff_n,Un_adh_t= glacier.F_result_typelaw1(Un,Udn,Fc[n],Un_adh,eps_Pen) 
        
        P_nalpha,Ff_nalpha,Un_adh_t = glacier.F_result_typelaw1(un_i,udn_i,Fc[n+1],Un_adh,eps_Pen)
        
        r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
        delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
        
        un_i1 = un_i + uiai*delta_uddn_i
        udn_i1 = udn_i + udiai*delta_uddn_i
        uddn_i1 = uddn_i + delta_uddn_i
        
        i = 0
        
        while i<100 and np.max(np.abs(udn_i1 - udn_i))>eps_HHT :
            
            un_i = un_i1
            udn_i = udn_i1
            uddn_i = uddn_i1
            
            Cn_i = np.zeros((glacier.Np,glacier.Np))
            for k in range(glacier.lim_ind_weert+1,glacier.Np):
                if np.abs(udn_i[k])>=glacier.ef :
                    Cn_i[k][k] = glacier.m*glacier.Cw*glacier.dl_l[k]*(np.abs(udn_i[k])**(glacier.m - 1))
                else :
                    Cn_i[k][k] = glacier.Cl*glacier.dl_l[k]
            
            Ms_i = M + mc*Cn_i + mk*K
        
            P_nalpha,Ff_nalpha,Un_adh = glacier.F_result_typelaw1(un_i,udn_i,Fc[n+1],Un_adh,eps_Pen)
            
            r_i = p_alpha*(P_n) + p_alpha1*(P_nalpha) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
            
            i +=1
        
        en = np.max(np.abs(udn_i1 - udn_i)) #residual error
        
        # state vector at n+1
        Un = un_i1
        Udn = udn_i1
        Uddn = uddn_i1
        
        # Un_adh = Un_adh_t
        
        # enregistrement des données pour affichage
        if n%Nt_r==0:
            if n%1000 == 0:
                print(n) 
            T_plot.append(T[n])
            Ut_plot.append( Un )
            Utd_plot.append( Udn )
            Utdd_plot.append( Uddn )
            
            Ut_adh_plot.append( Un_adh )
            Ftf_plot.append( np.sum(P_nalpha) )
            Ftsismique_plot.append(  np.sum(Ff_nalpha) )
            Ftsismique_map.append( Ff_nalpha )
            Niter_implicite.append(i)
            ErrResi.append(en)
        
    return(np.array(T_plot),np.array(Ut_plot), np.array(Utd_plot), np.array(Utdd_plot), np.array(Ut_adh_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map),np.array(Niter_implicite),np.array(ErrResi))

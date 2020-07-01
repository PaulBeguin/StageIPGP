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


## construction class Glacier pour initialisation des valeurs et fonctions
class Glacier:
    
# Paramètre de construction physique:
#   masse volumique glace : rho_glace [kg.m^-3]
#   masse volumique eau : rho_eau [kg.m^-3]
#   épaisseur du glacier : H [m]
#   pente du glacier : alpha_glacier [deg]
#   accélaration de la pesanteur : g [m.s^-2]
#   coeff de friction Weertman : Cw [Pa]
#   coeff de friction Budd : Cb [Pa]
#   coeff de friction Schoof : Cs [Pa]
#   exposant de vitesse de glissement : m
#   coeff de glissement max Schoof : Cmax [Pa]
#   loi utilisée : type_law =0 Weertman, =1 Budd, =2 Schoof, =3 Coulomb

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
<<<<<<< HEAD
    def __init__(self,rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ef,Np,Ltot,h_im):
=======
    def __init__(self,rho_glace,rho_eau,H,alpha_glacier,g,E,Cw,m,ef,Np,Ltot,h_im,Np_r):
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
        
        ## Paramètre de construction physique:
        # paramètres géométrique et physique élémentaires
        self.rho_glace=rho_glace #masse volumique glace
        self.rho_eau=rho_eau #masse volumique eau
        self.alpha_glacier=alpha_glacier #pente du glacier
        self.g=g #accélération pesanteur
        self.E=E #module d'young de la glace
        
        #definition de la loi de frottement
        self.Cw=Cw
        self.m=m
        self.ef=ef
        self.Cl=self.Cw*((self.ef)**(self.m-1))
        
        ## Construction élémentaire de la géométrie du glacier et discrétisation
        # discrétisation du glacier
        self.H=H #épaisseur du glacier
        self.Ltot = Ltot #longueur du glacier
        self.Np = Np #nombre de blocs
        
        ## Setting of the glacier geometrie
        self.l_def()
        
        ## Mass matrix and mass of the blocs array
        self.M_bloc_def()
        self.M_def()
        
        # définition de l'immergement et pression effective
        self.h_im = h_im
        self.h_im_l = np.cos(self.alpha_glacier) * h_im * np.ones(self.Np+1) - np.sin(self.alpha_glacier) * self.l #profondeur d'immergement au frontière entre les blocs
        
<<<<<<< HEAD
=======
        ## Initialisation des paramètres d'optimisation géométrique
        self.Np_r = Np_r
        self.lp_plot = np.arange(0,self.Np,self.Np_r) #liste d'indice blocs
        self.Np_plot = len(self.lp_plot) #longueur réduite de bloc de liste
        self.L_plot = np.linspace(0,self.Ltot,self.Np_plot) #position des blocs à afficher
        
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
        ## Construction of matrxi stiffness
        self.K_def() #enregistrement utile pour calcul de Fe
        
        ## calling for all the function 
        self.P_hydro_def()
        self.P_xpoids_def()
        self.P_zpoids_def()
        
        self.P_eff_def()
        
        self.Ud_eq_def()
        
        ## calcul of C damping linearized matrix at steady state
        self.C_lin_eq = self.C_lin_def(self.Ud_eq)
    
<<<<<<< HEAD
    ## stiffness matrix
=======
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
    def K_def(self):
        
        Ke = (self.E*self.H/self.dl) * np.array([[1,-1],[-1,1]])
        K_ = np.zeros((self.Np+1,self.Np+1))
        for i in range(self.Np):
            K_[i:i+2,i:i+2] = K_[i:i+2,i:i+2] + Ke
        
        self.K = K_
        
    ## pression hydrostatique sur les blocs
    def P_hydro_def(self):
        
        P_hydro_ = np.zeros(self.Np+1)
        for k in range(self.Np+1): # for regular element in the middle of the glacier
            h0 = self.h_im_l[k]
            h1 = self.h_im_l[k]
            if h0>0 :
                if h1>=0 :
                    P_hydro_[k] = self.rho_eau * self.g * self.dl_l[k] * (h0+h1)/2
                elif h1<0 :
                    dlk = self.dl_l[k] * h1/(h1-h0)
                    P_hydro_[k] = self.rho_eau * self.g * h1 * dlk/2
        
        self.P_hydro = P_hydro_
    
    ## pression poids composante normale au sol
    def P_zpoids_def(self):
        
        P_zpoids_ = -np.cos(self.alpha_glacier) * self.g * self.M_bloc
        
        self.P_zpoids = P_zpoids_
    
    ## pression poids composante tangentielle au sol - et valeur de frottement sur tout le glacier
    def P_xpoids_def(self):
        
        P_xpoids_ = -np.sin(self.alpha_glacier) * self.g * self.M_bloc
        self.F_stat = np.sum(P_xpoids_ )
        
        self.P_xpoids = P_xpoids_ 

    ## pression effective aux blocs 
    def P_eff_def(self):
        
        P_hydro_ = self.P_hydro
        P_zpoids_ = self.P_zpoids
        
        P_eff_=(np.abs(P_hydro_ + P_zpoids_) - (P_hydro_ + P_zpoids_))/2 #partie positive
        
        self.P_eff = P_eff_
    
    
    ## Initialisation of steady sliding state
    def Ud_eq_def(self):
        
        P_eff_ = self.P_eff
        P_xpoids_ = self.P_xpoids
        
        Ud_eq_ = np.zeros(self.Np+1)
        for k in range(self.Np+1):
            Ud_eq_[k] = (P_xpoids_[k]/(self.Cw*self.dl_l[k]))**(1/self.m)
        
        self.Ud_eq = Ud_eq_
    
    
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
            l_[k] = l_[k] + dl_l_[k]
        dl_l_[self.Np] = dl_/2
        l_[self.Np] = l_[self.Np-1] + dl_l_[self.Np]
        
        self.dl = dl_
        self.l = l_
        self.dl_l = dl_l_
    
    ## Mass for each bloc function
    def M_bloc_def(self):
        
        M_bloc_ = np.zeros(self.Np+1)
        for k in range(self.Np+1):
            M_bloc_[k] = self.H * self.rho_glace * self.dl_l[k]
        
        self.M_bloc = M_bloc_
    
    
    ## Mass matrix construction
    def M_def(self):
        
        M_ = np.diag(self.M_bloc)
        
        self.M = M_
    
    
    ## Linearised Matrix for damping
    def C_lin_def(self,Ud):
        
        C_v = np.zeros(self.Np+1)
        for k in range(self.Np+1):
            if np.abs(Ud[k])>=self.ef :
                C_v[k] = self.m*self.Cw*self.dl_l[k]*(np.abs(Ud[k])**(self.m - 1))
            else :
                C_v[k] = self.Cl*self.dl_l[k]
        
        C_ = np.diag(C_v)
        return (C_)
    
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
    
    ## Force élastique
    def Fe(self,U):
        
        # retourne la force élastique aux blocs
        # K=self.K , U=np.transpose(U), multiplication avec np.dot
        
        return(-np.dot(self.K,np.transpose(U)))
    
    ## Force de frottement
    def Ff(self,Ud):
        
        # retourne la force de frottement aux blocs 
        # Ud la vitesse, Ff vecteur force de frottement, self.Np longueur
        # self.type_law type de loi de frottement 
        # self.dl taille de bloc, self.Cw,self.Cb,self.Cs coefficient de friction, self.Cmax loi de Schoof, self.m exposant
        
        Ff_=np.zeros(self.Np+1)
        for k in range(self.Np+1):
            if np.abs(Ud[k])>self.ef :
                Ff_[k]=-np.sign(Ud[k])*self.Cw*self.dl_l[k]*(np.abs(Ud[k])**(self.m))
            else :
                Ff_[k]=-self.Cl*self.dl_l[k]*Ud[k]
        return(Ff_)
    
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

## Fonction pour intégration temporelle et simulation
<<<<<<< HEAD
def simu(glacier,Fc,dt,Nt,Nt_r):
=======
def simu(glacier,Fc,dt,Nt,li_plot):
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
    
    # définition des vecteurs de calcul et de l'état initial du système
    Ui = np.zeros(glacier.Np+1) # déplacement - nul
    Udi = glacier.Ud_eq # vitesse - prise à l'équilibre
    Fi = glacier.Fpertur(Fc[0]) # chargement - Fc[0]
    
    # définition du système optimisé - retourné et affiché
<<<<<<< HEAD
    Ut_plot = [[0 for k in range(glacier.Np+1)]] # displacement
    Udt_plot = [[Udi[k] for k in range(glacier.Np+1)]] # vitesse affichée
    Uddt_plot = [[0 for k in range(glacier.Np+1)]] # acceleration
    Ftf_plot = [0] # force
    Ftsismique_plot = [0] # force de frottement comparée au frottement statique
    Ftsismique_map = [[0 for k in range(glacier.Np+1)]]
=======
    Ut_plot = [] # déplacement affiché
    Udt_plot = [] # vitesse affichée
    Ftf_plot = [] # force de frottement affichée
    Ftsismique_plot = [] # force de frottement comparée au frottement statique
    Ftsismique_map = []
    Rm_plot = glacier.R_plot() # matrice de passage de Np dof à Np_plot déplacement affiché
    Rfm_plot = glacier.Rf_plot()
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
    
    # matrice du problème
    M = glacier.M #masse
    M_inv = np.linalg.inv(M) #calcul de l'inverse de masse
    
    # Boucle temporelle - Nt-1 itérations pour obtenir longueur finale de Ui,Fi et Udi Nt
    print(" Début d'intégration temporelle \n Rang enregistré: ") #message pour suivi du calcul
<<<<<<< HEAD
    for i in range(Nt-1):
        
        # intégration de Verlet
        Ui = Ui + dt*(Udi + (dt/2)*np.dot(M_inv,Fi))
        Fi1 = glacier.Fe(Ui) + glacier.Ff(Udi) + glacier.Fpertur(Fc[i]) + glacier.P_xpoids
        Udi = Udi + (dt/2)*np.dot(M_inv,(Fi + Fi1))

        Fi = Fi1
        
        # enregistrement des données pour affichage
        if i%Nt_r==0:
            print(i) 
        
        Ut_plot.append( Ui )
        Udt_plot.append( Udi )
        Uddt_plot.append( Uddi )
        Ftf_plot.append( np.sum(glacier.Ff(Udi)) )
        Ftsismique_plot.append( np.sum(glacier.Ff(Udi)) + glacier.F_stat )
        Ftsismique_map.append(glacier.Ff(Udi) + glacier.P_xpoids)
        
    return(np.array(Ut_plot), np.array(Udt_plot), np.array(Uddt_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map),np.array(Niter_implicite))
=======
    for i in range(Nt):
        
        # intégration de Verlet
        # déplacement à ti+1
        Ui = Ui + dt*(Udi + (dt/2)*np.dot(M_inv,Fi))
        # force à ti+1
        Fi1 = glacier.Fe(Ui) + glacier.Ff(Udi) + glacier.Fpertur(Fc[i]) + glacier.P_xpoids
        # vitesse à ti+1
        Udi = Udi + (dt/2)*np.dot(M_inv,(Fi + Fi1))
        # stockage force 
        Fi = Fi1
        
        # enregistrement des données pour affichage
        if i in li_plot :
            if i%10000==0 :
                print(i) #affichage de l'indice enregistré
            Upi = Ui - glacier.Ud_eq*(i+1)*dt
            Udpi = Udi - glacier.Ud_eq
            Ut_plot.append( np.dot(Rm_plot,Upi) )
            Udt_plot.append( np.dot(Rm_plot,Udpi) )
            Ftf_plot.append( np.sum(glacier.Ff(Udi)) )
            Ftsismique_plot.append( np.sum(glacier.Ff(Udi)) + glacier.F_stat )
            Ftsismique_map.append(np.dot(Rfm_plot,(glacier.Ff(Udi) + glacier.P_xpoids)))
        
    return(np.array(Ut_plot), np.array(Udt_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map))
>>>>>>> 891ea16fec04372fc3a038ee303b1d844175f6fe
    
## Fonction pour intégration temporelle et simulation
def simu_alpha(glacier,Fc,dt,Nt,alpha,eps,Nt_r):

    # définition des vecteurs de calcul et de l'état initial du système
    Un = np.zeros(glacier.Np+1) # displacement
    Udn = glacier.Ud_eq # vitesse 
    Uddn = np.zeros(glacier.Np+1) # acceleration
    
    # arrays to return and at plotting 
    Ut_plot = [[0 for k in range(glacier.Np+1)]] # displacement
    Udt_plot = [[Udn[k] for k in range(glacier.Np+1)]] # vitesse affichée
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
        
        Cn_i = glacier.C_lin_def(udn_i)
        Ms_i = M + mc*Cn_i + mk*K
        
        r_i = p_alpha*(glacier.Fpertur(Fc[n]) + glacier.P_xpoids + glacier.Ff(Udn) + glacier.Fe(Un)) + p_alpha1*(glacier.Fpertur(Fc[n+1]) + glacier.P_xpoids + glacier.Ff(udn_i) + glacier.Fe(un_i)) - np.dot(M,uddn_i)
        delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
        
        un_i1 = un_i + uiai*delta_uddn_i
        udn_i1 = udn_i + udiai*delta_uddn_i
        uddn_i1 = uddn_i + delta_uddn_i
        
        i = 0
        
        while i<100 and np.max(np.abs(udn_i1 - udn_i))>eps :
            
            un_i = un_i1
            udn_i = udn_i1
            uddn_i = uddn_i1
            
            Cn_i = glacier.C_lin_def(udn_i)
            Ms_i = M + mc*Cn_i + mk*K
            
            r_i = p_alpha*(glacier.Fpertur(Fc[n]) + glacier.P_xpoids + glacier.Ff(Udn) + glacier.Fe(Un)) + p_alpha1*(glacier.Fpertur(Fc[n+1]) + glacier.P_xpoids + glacier.Ff(udn_i) + glacier.Fe(un_i)) - np.dot(M,uddn_i)
            delta_uddn_i = np.dot(np.linalg.inv(Ms_i),r_i)
            
            un_i1 = un_i + uiai*delta_uddn_i
            udn_i1 = udn_i + udiai*delta_uddn_i
            uddn_i1 = uddn_i + delta_uddn_i
            
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
        Ftf_plot.append( np.sum(glacier.Ff(Udn)) )
        Ftsismique_plot.append( np.sum(glacier.Ff(Udn)) + glacier.F_stat )
        Ftsismique_map.append(glacier.Ff(Udn) + glacier.P_xpoids)
        
    return(np.array(Ut_plot), np.array(Udt_plot), np.array(Uddt_plot), np.array(Ftf_plot), np.array(Ftsismique_plot),np.array(Ftsismique_map),np.array(Niter_implicite))
    
## Fonction pour intégration temporelle et simulation
def simu_loop_dl(glacier,Fc,dt,Nt,li_plot,nu_plotk):
    
    # définition des vecteurs de calcul et de l'état initial du système
    Ui = np.zeros(glacier.Np+1) # déplacement - nul
    Udi = glacier.Ud_eq # vitesse - prise à l'équilibre
    Fi = glacier.Fpertur(Fc[0]) # chargement - Fc[0]
    
    # définition des listes pour maximum de déplacement, vitesse, force frottement et force sismique
    U_max = [0 for k in range(len(nu_plotk))]
    Ud_max = [0 for k in range(len(nu_plotk))]
    Ff_max = 0
    Fsismique_max = 0
    
    # matrice du problème
    M = glacier.M() #masse
    M_inv = np.linalg.inv(M) #calcul de l'inverse de masse
    
    # Boucle temporelle - Nt-1 itérations pour obtenir longueur finale de Ui,Fi et Udi Nt
    for i in range(Nt-1):
        
        # intégration de Verlet
        # déplacement à ti+1
        Ui = Ui + dt*(Udi + (dt/2)*np.dot(M_inv,Fi))
        # force à ti+1
        Fi1 = glacier.Fe(Ui) + glacier.Ff(Udi) + glacier.Fpertur(Fc[i]) + glacier.P_xpoids
        # vitesse à ti+1
        Udi = Udi + (dt/2)*np.dot(M_inv,(Fi + Fi1))
        # stockage force 
        Fi = Fi1
        
        # enregistrement des données pour affichage
        if i in li_plot :
            
            U_pertur = Ui-dt*(i+1)*glacier.Ud_eq
            Ud_pertur = Udi-glacier.Ud_eq
            Ff_result = np.sum(glacier.Ff(Udi))
            Fsismique_result = np.sum(glacier.Ff(Udi)) + glacier.F_stat
            
            # recherche du maximum - reprise des listes d'affichage avec même réduction
            for j in range(len(nu_plotk)):
                nu = nu_plotk[j]
                if U_max[j] < abs(U_pertur[nu]):
                    U_max[j] = abs(U_pertur[nu])
                if Ud_max[j] < abs(Ud_pertur[nu]):
                    Ud_max[j] = abs(Ud_pertur[nu])
            
            if Ff_max < Ff_result :
                Ff_max = Ff_result
            if Fsismique_max < Fsismique_result:
                Fsismique_max = Fsismique_result
        
    return(U_max,Ud_max,Ff_max,Fsismique_max)

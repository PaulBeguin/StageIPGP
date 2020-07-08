# -*- coding: utf-8 -*-
## Class for plotting the results

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp2d
import matplotlib.pyplot as mpl
import numpy as np
import os
import shutil
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] =10,7
matplotlib.rcParams.update({'font.size': 40})

class Plot_figures :
    
    def __init__(self,file_name_,work_path_,file_path_,save_path_,save_folder_,lu_plot_,Ttot_,Ltot_,Fc_):
        
        self.Ttot = Ttot_
        self.Ltot = Ltot_
        self.lu_plot = lu_plot_
        self.file_name = file_name_
        self.work_path = work_path_
        self.file_path = file_path_
        self.save_path = save_path_
        self.save_folder = save_folder_
        self.Fc = Fc_
      
        self.get_result()
        self.Np_plot = len(self.Ut[0])
        self.Nt_plot = len(self.Ut)
        self.T_plot = np.linspace(0,self.Ttot,self.Nt_plot)
        
        self.U_axis = [0,self.Ttot,-3.5,0.5]
        self.Up_axis = [0,self.Ttot,-0.5,1]
        self.Ud_axis = [0,self.Ttot,-60,20]
        self.F_axis = [0,self.Ttot,-50,50]
        
    
        self.get_Ut_max()
        
        # self.get_Fc_note()
        
        self.get_label()
        
        self.get_save_folder()
        
        # self.get_perturbation_Ut()
        mpl.rcParams['figure.figsize'] = 35,18
        
        # self.plot_displacement()
        # self.plot_speed()
        # self.plot_strains()
        # self.plot_Niter()
        
        # self.get_map()
        # self.get_perturbation_Ut()
        # self.map_displacement()
        # self.map_sismique()
    
    
    def get_result(self):
        
        file_result = np.load(os.path.join(self.file_path,self.file_name) )
        self.Ut = file_result['Ut']
        self.Utd = file_result['Utd']
        self.Utdd = file_result['Utdd']
        self.Niter_implicite = file_result['Niter_implicite']
        file_result.close()

    
    def get_Ut_max(self):
        
        Ut_max_ = 0
        Ut_ind_ = 0
        for k in range(len(self.Ut)):
            if Ut_max_<self.Ut[k][0] :
                Ut_max_ = self.Ut[k][0]
                Ut_ind_ = k
        
        self.T_Ut_max = self.Ttot*(Ut_ind_/len(self.Ut))
        self.Ut_max = Ut_max_
        
        print("Déplacement maximale : "+ "%.3f"%(self.Ut_max*1e+3) + " mm \n")
        print("A l'instant : "+ "%.1f"%(self.T_Ut_max) + " s \n")
    
    
    def get_label(self):
        
        nu_plot_=[0] 
        label_list_=['L= 0 m']
        for lu in self.lu_plot[1:len(self.lu_plot)-1]:
            nu_plot_.append(int(lu*self.Np_plot/self.Ltot))
            label_list_.append('L= '+ "%.0f"%lu +' m')
        
        nu_plot_.append(len(self.Ut[0])-1)
        label_list_.append('L= '+ "%.0f"%self.Ltot + ' m') 
        
        self.nu_plot = nu_plot_
        self.label_list = label_list_
        
        label_list_Ud_ = self.label_list
        self.label_list_Ud = label_list_Ud_
        
        self.label_list_F = ['Force de contact','Force de réaction']
    
    
    def get_save_folder(self):
        
        # os.chdir(self.save_path)
        # 
        # if self.save_folder in os.listdir() :
        #     shutil.rmtree(self.save_folder)
        #     
        # os.makedirs(self.save_folder)
        os.chdir(os.path.join(self.save_path,self.save_folder))
    
    
    def plot_displacement(self):
        
        name_fig1 = "Deplacement"
        fg1 = plt.figure(1)
        
        for i in range(len(self.nu_plot)):
            U_nu_plot = []
            nu = self.nu_plot[i]
            for k in range(self.Nt_plot):
                U_nu_plot.append(self.Ut[k][nu]*1000)
            plt.plot(self.T_plot,U_nu_plot)
        
        plt.legend(self.label_list,ncol=3)
        # plt.axis(self.U_axis)
        
        plt.title('Deplacement total')
        plt.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('$U_x$ (mm)')
        plt.draw()
        plt.tight_layout() 
        fg1.savefig(name_fig1+".svg")
    
    
    def plot_speed(self):
        
        name_fig2 = "Vitesse"
        fg2 = plt.figure(2)
        ax2 = fg2.gca()
        
        for i in range(len(self.nu_plot)):
            Ud_nu_plot = []
            nu = self.nu_plot[i]
            for k in range(self.Nt_plot):
                Ud_nu_plot.append(self.Utd[k][nu]*1e+6)
            ax2.plot(self.T_plot,Ud_nu_plot)
        

        
        ax2.legend(self.label_list_Ud,ncol=3)
        # ax2.axis(self.Ud_axis)
        
        plt.title('Vitesse de perturbation du glacier en fonction du temps')
        ax2.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('$\dot{U}_x$ ($\mu$ m.s$^{-1}$)')
        plt.draw()
        plt.tight_layout()
        fg2.savefig(name_fig2+'.svg')
    
    
    # def plot_strains(self):
    #     
    #     name_fig3 = "Fsismique"
    #     fg3 = plt.figure(3)
    #     ax3 = fg3.gca()
    #     
    #     ax3.plot(self.T_plot,self.Fc*1e-6)
    #     ax3.plot(self.T_plot[1:self.Nt_plot],self.Ftsismique[1:self.Nt_plot]*1e-6)

   ##       ax3.legend(self.label_list_F,ncol=2)
    #     ax3.axis(self.F_axis)
    #     
    #     plt.title("Force de contact entre glacier et sol comparée a celle de l'iceberg sur le glacier")
    #     ax3.grid(True)
    #     plt.xlabel('Temps $(s)$')
    #     plt.ylabel('Force $e_{x}$ $(MN.m^{-1})$')
    #     plt.draw()
    #     fg3.savefig(name_fig3+'.svg')
    
    
    def plot_Niter(self):
        
        name_fig4 = "NiterHHT"
        fg4 = plt.figure(4)
        ax4 = fg4.gca()
        
        ax4.plot(self.T_plot[:-1],self.Niter_implicite)
        
        plt.title("Nombre d'iteration a chaque pas de temps - schema HHT")
        plt.xlabel('Instant simule - temps $t$ (s)')
        plt.ylabel('$N_{iter}$')
        plt.grid(True)
        plt.draw()
        plt.tight_layout()
        fg4.savefig(name_fig4+'.svg')
    
    
    def get_map(self):
        
        self.L_map = np.linspace(0,self.Ltot,self.Np_plot)
        self.T_mesh , self.L_mesh = np.meshgrid(self.T_plot,self.L_map)
    
    
    def map_displacement(self):
        
        name_fig6 = "MapDisplacement"
        fg6 = plt.figure(6)
        ax6 = fg6.gca()
        
        im6 = ax6.pcolormesh(self.T_mesh,self.L_mesh,np.transpose(self.Ut)*1000,shading='gouraud')
        
        fg6.colorbar(im6 , ax=ax6)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('Position du bloc de glacier $x$ (km)')
        plt.title('Carte de la perturbation du deplacement $\Delta U_x$ (mm) dans le glacier en fonction du temps et de la position du bloc')
        #plt.show()
        plt.tight_layout()
        fg6.savefig(name_fig6+".png")
    
     
    def map_sismique(self):
         
        name_fig7 = "MapSismique"
        fg7 = plt.figure(7)
        ax7 = fg7.gca()
         
        im7 = ax7.pcolormesh(self.T_mesh,self.L_mesh*1e-3,np.transpose(self.Ftsismique_map)*1e-6,shading='gouraud')
        plt.colorbar(im7 , ax=ax7)
         
        plt.xlabel('Temps $(s)$')
        plt.ylabel('Position du bloc de glacier $(km)$')
        plt.title('Force sismique $(MN.m^{-1})$ dans le glacier en fonction du temps et de la position du bloc')
        #plt.show()
        fg7.savefig(name_fig7+".png")


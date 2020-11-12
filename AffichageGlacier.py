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

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] =35,18
matplotlib.rcParams.update({'font.size': 12})

class Plot_figures :
    
    def __init__(self,file_name_,file_Fcname_,work_path_,file_path_,save_path_,save_folder_,lu_plot_,Ttot_,Ltot_,Fc_,ud_reg_):
        
        self.Ttot = Ttot_
        self.Ltot = Ltot_
        self.lu_plot = lu_plot_
        self.file_name = file_name_
        self.file_Fcname = file_Fcname_
        self.work_path = work_path_
        self.file_path = file_path_
        self.save_path = save_path_
        self.save_folder = save_folder_
        self.Fc = Fc_
        self.ud_reg = ud_reg_
        
        self.get_result()
        self.Np_plot = len(self.Ut[0])
        self.Nt_plot = len(self.Ut)
        self.T_plot = np.linspace(0,self.Ttot,self.Nt_plot)
        
        self.U_axis = [0,self.Ttot,-3.5,0.5]
        self.Up_axis = [0,self.Ttot,-0.5,1]
        self.Ud_axis = [0,self.Ttot,-60,20]
        self.F_axis = [0,self.Ttot,-50,50]
        
        self.title_font = {'size':'20'} 
        self.label_font = {'size':'20'}
        self.legend_font = {'size':'15'}
        
        self.get_Ut_max()
        
        self.get_Fc_note()
        self.T_retour_plot = np.array([self.T_retour,self.T_retour])
        self.T_max_plot = np.array([self.T_Fcmax,self.T_Fcmax])
        
        self.U_retour_plot = self.U_axis[2:4]
        self.Up_retour_plot = self.Up_axis[2:4]
        self.Ud_retour_plot = self.Ud_axis[2:4]
        self.F_retour_plot = self.F_axis[2:4]
        
        self.U_max_plot = self.U_retour_plot
        self.Up_max_plot = self.Up_retour_plot
        self.Ud_max_plot = self.Ud_retour_plot
        self.F_max_plot = self.F_retour_plot
        
        self.Ud_eq_plot = np.array([self.Utd[0][0],self.Utd[0][0]])
        self.Ud_reg_plot = np.array([self.ud_reg,self.ud_reg])
        self.Ud_reg_plot2 = np.array([-self.ud_reg,-self.ud_reg])
        self.T_eq_plot = self.Ud_axis[0:2]
        self.T_reg_plot = self.T_eq_plot
        
        self.get_label()
        
        self.get_save_folder()
        
        self.get_perturbation_Ut()
        # mpl.rcParams['figure.figsize'] = 35,18
        
        # self.plot_displacement()
        # self.plot_speed()
        # self.plot_strains()
        # self.plot_Niter()
        
        self.get_map()
        self.get_perturbation_Ut()
        # self.map_displacement()
        # self.map_sismique()
    
    
    def get_result(self):
        
        file_result = np.load(self.file_path + "\\" + self.file_name )
        self.Ut = file_result['Ut']
        self.Utd = file_result['Utd']
        self.Utdd = file_result['Utdd']
        self.Ftf = file_result['Ftf']
        self.Ftsismique = file_result['Ftsismique']
        self.Ftsismique_map = file_result['Ftsismique_map']
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
    
    
    def get_Fc_note(self):
        
        file_Fc = open(self.work_path + '\\' + self.file_Fcname,'r')
        file_Fc_txt = file_Fc.read()
        file_Fc.close()
        file_Fc_txt_lines = file_Fc_txt.splitlines()
        
        self.Nt_retour = int(file_Fc_txt_lines[0])
        self.T_retour = float(file_Fc_txt_lines[1])
        self.Fc_max = float(file_Fc_txt_lines[2])
        self.T_Fcmax = float(file_Fc_txt_lines[3])
    
    
    def get_label(self):
        
        nu_plot_=[0] 
        label_list_=['L= 0 m']
        for lu in self.lu_plot[1:len(self.lu_plot)-1]:
            nu_plot_.append(int(lu*self.Np_plot/self.Ltot))
            label_list_.append('L= '+ "%.0f"%lu +' m')
        
        nu_plot_.append(len(self.Ut[0])-1)
        label_list_.append('L= '+ "%.0f"%self.Ltot + ' m') 
        label_list_.append('Contact maximal')
        label_list_.append('Retournement complet')
        
        self.nu_plot = nu_plot_
        self.label_list = label_list_
        
        label_list_Ud_ = label_list_.copy()
        label_list_Ud_.append("Vitesse d'équilibre")
        label_list_Ud_.append('Limite de régularisation')
        self.label_list_Ud = label_list_Ud_
        
        self.label_list_F = ['Force de contact','Force sismique','Contact maximal','Retournement complet']
    
    
    def get_save_folder(self):
        
        # os.chdir(self.save_path)
        # 
        # if self.save_folder in os.listdir() :
        #     shutil.rmtree(self.save_folder)
        #     
        # os.makedirs(self.save_folder)
        os.chdir(self.save_path + '\\' + self.save_folder)
    
    
    def plot_displacement(self):
        
        name_fig1 = "Deplacement"
        fg1 = plt.figure(1)
        
        plt.subplot(121)
        
        for i in range(len(self.nu_plot)):
            U_nu_plot = []
            nu = self.nu_plot[i]
            for k in range(self.Nt_plot):
                U_nu_plot.append(self.Ut[k][nu]*1000)
            plt.plot(self.T_plot,U_nu_plot)
        
        plt.plot(self.T_max_plot , self.U_max_plot , 'k--',linewidth=3)
        plt.plot(self.T_retour_plot , self.U_retour_plot , 'k--',linewidth=1)
        plt.legend(self.label_list,ncol=3)
        plt.axis(self.U_axis)
        
        plt.title('Déplacement total',self.title_font)
        plt.grid(True)
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('$U_x$ (mm)',self.label_font)

        plt.subplot(122)
        
        for i in range(len(self.nu_plot)):
            Up_nu_plot = []
            nu = self.nu_plot[i]
            for k in range(self.Nt_plot):
                Up_nu_plot.append(self.Utp[k][nu]*1000)
            plt.plot(self.T_plot,Up_nu_plot)
        
        plt.plot(self.T_max_plot , self.Up_max_plot , 'k--',linewidth=3)
        plt.plot(self.T_retour_plot , self.Up_retour_plot , 'k--',linewidth=1)
        plt.legend(self.label_list,ncol=3)
        plt.axis(self.Up_axis)
        
        plt.title('Déplacement perturbé',self.title_font)
        plt.grid(True)
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('$U_x$ (mm)',self.label_font)
        
        # plt.show()
        

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
        
        ax2.plot(self.T_max_plot , self.Ud_max_plot , 'k--',linewidth=3)  
        ax2.plot(self.T_retour_plot , self.Ud_retour_plot , 'k--',linewidth=1)
        ax2.plot(self.T_eq_plot , self.Ud_eq_plot*1e+6 , 'g--', linewidth=2)
        ax2.plot(self.T_reg_plot , self.Ud_reg_plot*1e+6 , 'r--',linewidth=1)
        ax2.plot(self.T_reg_plot , self.Ud_reg_plot2*1e+6 , 'r--',linewidth=1)
        
        ax2.legend(self.label_list_Ud,ncol=3)
        ax2.axis(self.Ud_axis)
        
        plt.title('Vitesse de perturbation du glacier en fonction du temps',self.title_font)
        ax2.grid(True)
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('$\dot{U}_x$ ($\mu$ m.s$^{-1}$)',self.label_font)
        plt.draw()
        fg2.savefig(name_fig2+'.svg')
        # fg2.savefig(name_fig2+'.svg')
    
    def plot_strains(self):
        
        name_fig3 = "Fsismique"
        fg3 = plt.figure(3)
        ax3 = fg3.gca()
        
        ax3.plot(self.T_plot,self.Fc*1e-6)
        ax3.plot(self.T_plot[1:self.Nt_plot],self.Ftsismique[1:self.Nt_plot]*1e-6)

        ax3.plot(self.T_max_plot , self.F_max_plot , 'k--',linewidth=3)
        ax3.plot(self.T_retour_plot , self.F_retour_plot , 'k--',linewidth=1)
        ax3.legend(self.label_list_F,ncol=2)
        ax3.axis(self.F_axis)
        
        plt.title("Force de contact entre glacier et sol comparee a celle de l'iceberg sur le glacier",self.title_font)
        ax3.grid(True)
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('Force $F_x$ (MN.m$^{-1}$)',self.label_font)
        plt.draw()
        fg3.savefig(name_fig3+'.svg')
        # fg3.savefig(name_fig3+'.svg')

    
    def plot_Niter(self):
        
        name_fig4 = "NiterHHT"
        fg4 = plt.figure(4)
        ax4 = fg4.gca()
        
        ax4.plot(self.T_plot[:-1],self.Niter_implicite)
        
        plt.title("Nombre d'iteration a chaque pas de temps - schema HHT",self.title_font)
        plt.xlabel('Instant simule - temps $t$ (s)',self.label_font)
        plt.ylabel('$N_{iter}$',self.label_font)
        plt.grid(True)
        plt.draw()
        fg4.savefig(name_fig4+'.svg')
        # fg4.savefig(name_fig4+'.svg')
    
    def get_map(self):
        
        self.L_map = np.linspace(0,self.Ltot,self.Np_plot)
        self.T_mesh , self.L_mesh = np.meshgrid(self.T_plot,self.L_map)
    
    def get_perturbation_Ut(self):
        
        Ut_eq_ = np.zeros((self.Nt_plot,self.Np_plot))
        for k in range(self.Nt_plot):
            for i in range(self.Np_plot):
                Ut_eq_[k][i] = self.T_plot[k]*self.Utd[0][i]
        
        Utp_ = self.Ut - Ut_eq_
        self.Utp = Utp_
    
    def map_displacement(self):
        
        name_fig6 = "MapDisplacement"
        fg6 = plt.figure(6)
        ax6 = fg6.gca()
        
        # levels = MaxNLocator(nbins=15).tick_values(self.Utp.min(), self.Utp.max())
        # 
        # cmap = plt.get_cmap('PiYG')
        # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        
        im6 = ax6.pcolormesh(self.T_mesh,self.L_mesh,np.transpose(self.Utp)*1000,shading='gouraud')
        
        fg6.colorbar(im6 , ax=ax6)
        plt.xlabel('Temps $r$ (s)',self.label_font)
        plt.ylabel('Position du bloc de glacier $x$ (km)',self.label_font)
        plt.title('Carte de la perturbation du déplacement(mm) dans le glacier en fonction du temps et de la position du bloc',self.title_font)
        plt.show()
        fg6.savefig(name_fig6+".svg")
        # fg6.savefig(name_fig6+".svg")
    
    def map_displacement2(self):
        
        name_fig8 = "MapDisplacement2"
        fg8 = plt.figure(8)
        ax8 = fg8.gca()
        
        f = interp2d(self.T_mesh,self.L_mesh,np.transpose(self.Utp)*1000,kind='cubic')
        self.T_mesh2 = np.linspace(0,self.Ttot,self.Nt_plot)
        self.L_mesh2 = np.linspace(0,self.Ltot,self.Np_plot)
        
        data1 = f(self.T_mesh2,self.L_mesh2)
        X2, Y2 = np.meshgrid(self.T_mesh2,self.L_mesh2)

        im8 = ax8.pcolormesh(X2, Y2, data1, cmap='RdBu')
        fg8.colorbar(im8 , ax=ax8)
        
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('Position du bloc de glacier $x$ (km)',self.label_font)
        plt.title('Carte de la perturbation du deplacement(mm) dans le glacier en fonction du temps et de la position du bloc',self.title_font)
        plt.show()
        fg8.savefig(name_fig8+".svg")
        # fg8.savefig(name_fig8+".svg")
    
    def map_sismique(self):
        
        name_fig7 = "MapSismique"
        fg7 = plt.figure(7)
        ax7 = fg7.gca()
        
        im7 = ax7.pcolormesh(self.T_mesh,self.L_mesh*1e-3,np.transpose(self.Ftsismique_map)*1e-6,shading='gouraud')
        plt.colorbar(im7 , ax=ax7)
        
        plt.xlabel('Temps $t$ (s)',self.label_font)
        plt.ylabel('Position du bloc de glacier $x$ (km)',self.label_font)
        plt.title('Force sismique (MN.m^${-1}$) dans le glacier en fonction du temps et de la position du bloc',self.title_font)
        plt.show()
        fg7.savefig(name_fig7+".svg")
        # fg7.savefig(name_fig7+".svg")
    

class Plot_figure_dl_loop :
    
    def __init__(self,file_Fcname_,work_path_,results_path_,results_file_,save_path_,save_folder_,dl_list_,Ttot_,ud_reg_):
        
        self.file_Fcname = file_Fcname_
        self.work_path = work_path_
        self.results_path = results_path_
        self.results_file = results_file_
        self.save_path = save_path_
        self.save_folder = save_folder_
        self.dl_list = dl_list_
        self.Ttot = Ttot_
        self.ud_reg = ud_reg_
        
        self.get_result()
        self.N_dl = len(self.U)
        self.Nt_plot = len(self.U[0])
        self.T_plot = np.linspace(0,self.Ttot,self.Nt_plot)
        self.Up = self.T_plot*self.Ud[0][0]
        
        self.U_axis = [0,self.Ttot,-3.5,0.5]
        self.Up_axis = [0,self.Ttot,-0.5,1]
        self.Ud_axis = [0,self.Ttot,-60,20]
        self.F_axis = [0,self.Ttot,-150,50] 
        
        self.U_error_axis = [1000,10,0,0.025]
        self.Ud_error_axis = [1000,10,0,0.0007]
        self.Fsis_error_axis = [1000,10,0,2.5*1e+7]
        
        self.get_Fc_note()
        self.T_retour_plot = np.array([self.T_retour,self.T_retour])
        self.T_max_plot = np.array([self.T_Fcmax,self.T_Fcmax])
        
        self.U_retour_plot = self.U_axis[2:4]
        self.Up_retour_plot = self.Up_axis[2:4]
        self.Ud_retour_plot = self.Ud_axis[2:4]
        self.F_retour_plot = self.F_axis[2:4]
        
        self.U_max_plot = self.U_retour_plot
        self.Up_max_plot = self.Up_retour_plot
        self.Ud_max_plot = self.Ud_retour_plot
        self.F_max_plot = self.F_retour_plot
        
        self.Ud_reg_plot = np.array([self.ud_reg,self.ud_reg])
        self.Ud_reg_plot2 = np.array([-self.ud_reg,-self.ud_reg])
        self.T_reg_plot = self.Ud_axis[0:2]
        
        self.get_label()
        
        self.get_save_folder()
        
        self.get_error()
        
        # self.plot_displacement()
        # self.plot_speed()
    
    def get_result(self):
        
        file_result = np.load(self.results_path + "\\" + self.results_file )
        self.U = file_result['U']
        self.Ud = file_result['Ud']
        self.Fsis = file_result['Fsis']
        file_result.close()
    
    
    def get_Fc_note(self):
        
        file_Fc = open(self.work_path + '\\' + self.file_Fcname,'r')
        file_Fc_txt = file_Fc.read()
        file_Fc.close()
        file_Fc_txt_lines = file_Fc_txt.splitlines()
        
        self.Nt_retour = int(file_Fc_txt_lines[0])
        self.T_retour = float(file_Fc_txt_lines[1])
        self.Fc_max = float(file_Fc_txt_lines[2])
        self.T_Fcmax = float(file_Fc_txt_lines[3])
    
    
    def get_label(self):
        
        label_list_ = []
        for dl in self.dl_list:
            label_list_.append("dl = " + str(dl))
        
        label_list_.append('Contact maximal')
        label_list_.append('Retournement complet')
        self.label_list_U = label_list_
        
        label_list_Ud_ = label_list_
        label_list_Ud_.append('Limite de régularisation')
        self.label_list_Ud = label_list_Ud_
        
        self.label_list_F = self.label_list_U
    
    
    def get_save_folder(self):
        
        # os.chdir(self.save_path)
        # 
        # if self.save_folder in os.listdir() :
        #     shutil.rmtree(self.save_folder)
        #     
        # os.makedirs(self.save_folder)
        os.chdir(self.save_path + '\\' + self.save_folder)
    
    
    def plot_displacement(self):
        
        name_fig1 = "Deplacement"
        fg1 = plt.figure(1)
        
        plt.subplot(121)
        
        for i in range(self.N_dl):
            U_plot = []
            for k in range(self.Nt_plot):
                U_plot.append(self.U[i][k]*1000)
            plt.plot(self.T_plot,U_plot)
        
        plt.plot(self.T_max_plot , self.U_max_plot , 'k--',linewidth=3)
        plt.plot(self.T_retour_plot , self.U_retour_plot , 'k--',linewidth=1)
        plt.legend(self.label_list_U,ncol=3)
        plt.axis(self.U_axis)
        
        plt.title('Deplacement total')
        plt.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('$U_x$ (mm)')
        plt.draw()
        
        plt.subplot(122)
        
        for i in range(len(self.dl_list)):
            Up_plot = []
            for k in range(self.Nt_plot):
                Up_plot.append((self.U[i][k]-self.Up[k])*1000)
            plt.plot(self.T_plot,Up_plot)
        
        plt.plot(self.T_max_plot , self.Up_max_plot , 'k--',linewidth=3)
        plt.plot(self.T_retour_plot , self.Up_retour_plot , 'k--',linewidth=1)
        plt.legend(self.label_list_U,ncol=3)
        plt.axis(self.Up_axis)
        
        plt.title('Deplacement perturbe')
        plt.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('$U_x$ (mm)')
        plt.draw()
        
        fg1.savefig(name_fig1, bbox_inches=None)
    
    
    def plot_speed(self):
        
        name_fig2 = "Vitesse"
        fg2 = plt.figure(2)
        ax2 = fg2.gca()
        
        for i in range(len(self.dl_list)):
            Ud_plot = []
            for k in range(self.Nt_plot):
                Ud_plot.append(self.Ud[i][k]*1e+6)
            ax2.plot(self.T_plot,Ud_plot)
        
        ax2.plot(self.T_max_plot , self.Ud_max_plot , 'k--',linewidth=3)  
        ax2.plot(self.T_retour_plot , self.Ud_retour_plot , 'k--',linewidth=1)
        ax2.plot(self.T_reg_plot , self.Ud_reg_plot*1e+6 , 'r--',linewidth=1)
        ax2.plot(self.T_reg_plot , self.Ud_reg_plot2*1e+6 , 'r--',linewidth=1)
        
        ax2.legend(self.label_list_Ud,ncol=3)
        ax2.axis(self.Ud_axis)
        
        plt.title('Vitesse de perturbation du glacier en fonction du temps')
        ax2.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('$\dot{U}_x$ ($\mu$ m.s$^{-1}$)')
        plt.draw()
        fg2.savefig(name_fig2, bbox_inches=None)
    
    
    def plot_strains(self):
        
        name_fig3 = "Fsismique"
        fg3 = plt.figure(3)
        ax3 = fg3.gca()
        
        for i in range(self.N_dl):
            Fsis_i = []
            for k in range(self.Nt_plot):
                Fsis_i.append(self.Fsis[i][k]/self.dl_list[i]*1e-3)
            ax3.plot(self.T_plot,Fsis_i)

        ax3.plot(self.T_max_plot , self.F_max_plot , 'k--',linewidth=3)
        ax3.plot(self.T_retour_plot , self.F_retour_plot , 'k--',linewidth=1)
        ax3.legend(self.label_list_F,ncol=2)
        ax3.axis(self.F_axis)
        
        plt.title("Force de contact entre glacier et sol par longueur de bloc")
        ax3.grid(True)
        plt.xlabel('Temps $t$ (s)')
        plt.ylabel('Force $F_{x}$ (kN.m$^{-2}$)')
        plt.draw()
        fg3.savefig(name_fig3,bbox_inches=None)
    
    
    def get_error(self):
        
        U_error_, Ud_error_, Fsis_error_ = [], [], []
        
        for i in range(self.N_dl-1):
            
            U_error_i_ = np.trapz(np.abs(self.U[-1]-self.U[i]),self.T_plot)/self.Ttot
            Ud_error_i_ = np.trapz(np.abs(self.Ud[-1]-self.Ud[i]),self.T_plot)/self.Ttot
            Fsis_error_i_ = np.trapz(np.abs((self.Fsis[-1]/self.dl_list[-1])-(self.Fsis[i]/dl_list[i])),self.T_plot)/self.Ttot
            
            U_error_.append(U_error_i_)
            Ud_error_.append(Ud_error_i_)
            Fsis_error_.append(Fsis_error_i_)
        
        self.U_error = U_error_
        self.Ud_error = Ud_error_
        self.Fsis_error = Fsis_error_
    
    
    def plot_error(self):
        
        name_fig4 = "ErrorDlLoop"
        fg4 = plt.figure(4)
        ax4 = fg4.gca()
        
        plt.subplot(131)
        plt.semilogx()
        plt.plot(self.dl_list[0:self.N_dl-1],self.U_error,'o--')
        plt.axis(self.U_error_axis)
        plt.title('Erreur en deplacement')
        plt.grid(True)
        plt.xlabel('$\delta_l$ (m)')
        plt.ylabel('$U_{x_{error}}$ (m)')
        plt.draw()
        
        plt.subplot(132)
        plt.semilogx()
        plt.plot(self.dl_list[0:self.N_dl-1],self.Ud_error,'o--')
        plt.axis(self.Ud_error_axis)
        plt.title('Erreur en vitesse')
        plt.grid(True)
        plt.xlabel('$\delta_l$ (m)')
        plt.ylabel('$Ud_{x_{error}}$ (m.s$^{-1}$)')
        plt.draw()
        
        plt.subplot(133)
        plt.semilogx()
        plt.plot(self.dl_list[0:self.N_dl-1],self.Fsis_error,'o--')
        plt.axis(self.Fsis_error_axis)
        plt.title('Erreur en vitesse')
        plt.grid(True)
        plt.xlabel('$\delta_l$ (m)')
        plt.ylabel('$Fsis_{x_{error}}$ (N.m$^{-2}$)')
        plt.draw()
        fg4.savefig(name_fig4,bbox_inches=None)

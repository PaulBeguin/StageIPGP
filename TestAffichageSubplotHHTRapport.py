# -*- coding: utf-8 -*-
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

resultsfile0 = "ResultsTestHHT1bloc.npz"
resultsfile1 = "ResultsTestHHT1bloc2.npz"
resultsfile2 = "ResultsTestHHT1bloc3.npz"
resultsfile3 = "ResultsTestHHT1bloc4.npz"

file_results_reloaded0 = np.load(path_resultsfile + "\\" + resultsfile0)
U0 = file_results_reloaded0['U']
Ud0 = file_results_reloaded0['Ud']
F0 = file_results_reloaded0['F']
file_results_reloaded0.close()

file_results_reloaded1 = np.load(path_resultsfile + "\\" + resultsfile1)
U1 = file_results_reloaded1['U']
Ud1 = file_results_reloaded1['Ud']
F1 = file_results_reloaded1['F']
file_results_reloaded1.close()

file_results_reloaded2 = np.load(path_resultsfile + "\\" + resultsfile2)
U2 = file_results_reloaded2['U']
Ud2 = file_results_reloaded2['Ud']
F2 = file_results_reloaded2['F']
file_results_reloaded2.close()

file_results_reloaded3 = np.load(path_resultsfile + "\\" + resultsfile3)
U3 = file_results_reloaded3['U']
Ud3 = file_results_reloaded3['Ud']
F3 = file_results_reloaded3['F']
file_results_reloaded3.close()

plt.figure(1)
plt.subplot(131)
T_plot = np.linspace(0,Ttot,len(U))
plt.plot(T_plot0,U0*1e+3)
plt.plot(T_plot3,U3*1e+3)
plt.plot(T_plot3,ud_eq*T_plot3*1e+3,'g--')
plt.legend(['Deplacement dt_reg','Deplacement dt_reg*1000','Equilibre statique'])
plt.grid(True)
plt.xlabel('Temps (s)')
plt.ylabel('Deplacement ($mm $)')
plt.gca().set_title('Déplacement du bloc')

plt.subplot(132)
plt.plot(T_plot0,Ud0*1e+6)
plt.plot(T_plot3,Ud3*1e+6)
plt.plot(T_plot3,ud_eq*np.ones((len(T_plot3)))*1e+6,'g--')
plt.legend(['Vitesse dt_reg','Vitesse dt_reg*1000','Equilibre'])
plt.grid(True)
plt.xlabel('Temps (s)')
plt.ylabel('Vitesse ($\mu m .s^{-1}$)')
plt.gca().set_title('Vitesse du bloc')

plt.subplot(133)
plt.plot(T_plot0,F0)
plt.plot(T_plot3,F3)
plt.legend(['Force dt_reg','Force dt_reg*1000'])
plt.grid(True)
plt.xlabel('Temps (s)')
plt.ylabel('Forces (N)')
plt.gca().set_title('Force sismique')

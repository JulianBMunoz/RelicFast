######################################################
####         USER INPUTS
######################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()

######################################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"



omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

m_ax = np.array([
    10**-25, 
    10**-26, 
    10**-27, 
    10**-28, 
    10**-29, 
    10**-30, 
    np.power(10., -31),
    np.power(10., -32)
])
omega_ax = np.array([
    1.0e-12*omega_cdm_LCDM, 
    0.01*omega_cdm_LCDM, 
    0.02*omega_cdm_LCDM,
    0.03*omega_cdm_LCDM,
    0.04*omega_cdm_LCDM,
    0.05*omega_cdm_LCDM,
    0.06*omega_cdm_LCDM,
    0.07*omega_cdm_LCDM,
    0.08*omega_cdm_LCDM,
    0.09*omega_cdm_LCDM,
    0.10*omega_cdm_LCDM
])

b1l_kmax = np.zeros((len(m_ax), len(omega_ax)))
b1e_kmax = np.zeros((len(m_ax), len(omega_ax)))

b1l_max = np.zeros((len(m_ax), len(omega_ax)))
b1e_max = np.zeros((len(m_ax), len(omega_ax)))

sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-4
kmax = 0.5
Nk = 50

######################################################



######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

#data_pm = []
#data_tf = []
#data_eulbias = []
#data_lagbias = []
#axion_background = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
            "#define boltzmann_tag  _CLASS_", "#define boltzmann_tag _AXIONCAMB_"
    ).replace(
        "#define boltzmann_tag  _CAMB_", "#define boltzmann_tag _AXIONCAMB_"
    )
        
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

os.chdir(rfpath)

# Run RelicFast for each axion abundance 
for ax_idx, ax_val in enumerate(m_ax): 
    for o_idx, o_val in enumerate(omega_ax):
        print("Running RelicFast + axionCAMB for m_ax = "+f"{ax_val:.2e}"+", omega_ax = "+f"{o_val:.3f}")
        
        omega_nu = sum_massive_nu/93.2
        Mnu = sum_massive_nu # Units: eV
        
        f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
        f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
        
        omega_nu = Mnu/93.14
        
        omega_cdm = omega_cdm_LCDM - o_val
        f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
        f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)
        
        h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+o_val)/(omega_b_LCDM+omega_cdm_LCDM))
        #h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
        
        kfs = np.pi * np.sqrt(ax_val*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)
        
        # Clear old data 
        os.system('rm -r '+rfpath_outputsuffix)
        os.system('rm -r '+'Boltzmann_0')
        os.system('rm -r '+'Boltzmann_1')
        os.system('rm -r '+'Boltzmann_2')
        os.system('rm ./run.ini')
        
        # RUN RelicFast + solver for each neutrino mass choice 
        reading_file = open("2006.09395.ini", "r")
        new_file_content = ""
        for line in reading_file:
            stripped_line = line.strip()
            new_line = stripped_line.replace(
                "mnu1 = 0.0", "mnu1 = "+f'{Mnu:.2e}'
            ).replace(
                "mnu2 = 0.0", "mnu2 = "+f'{Mnu:.2e}'
            ).replace(
                "m_SN = 0.02", "m_SN = "+f'{Mnu:.2e}'
            ).replace(
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{o_val:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm:.6e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{ax_val:.6e}'
            )
            
            
            if (Mnu!=0.): 
                new_line = new_line.replace(
                    "tag_sterile_nu = 0", "tag_sterile_nu = 1"
                )
                
            new_file_content += "    " + new_line +"\n"
        reading_file.close()
        writing_file = open("run.ini", "w")
        writing_file.write(new_file_content)
        writing_file.close()
        
        os.system('./relicfast run.ini')
        
        # Collect axion background evolution
        #axion_background.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_background.dat"))    
        
        # Collect power spectrum for requested redshift 
        output_dirs = os.listdir(rfpath_boltzmannsuffix)
        output_dirs.remove('_params.ini')
        
        #z_vals = np.array([])
        #
        #for str_idx, str_val in enumerate(output_dirs):
        #    if (str_val[0:15]=='_transfer_out_z'): 
        #        z_vals = np.append(
        #            z_vals, 
        #            np.float(str_val.split('_transfer_out_z')[1])
        #        )
        #            
        #z_vals = np.sort(z_vals)
            
    
        #pm_idx = np.argmin(np.abs(z_vals - redshift))
        #pm_idx = np.argmin(np.abs(z_vals - 0.0))
        #print('Requested/found redshift: ', redshift, z_vals[pm_idx])
        
        #print('Loading: ')
        
        #print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        #print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        
        #data_pm.append(
        #    np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        #)
        
        #data_tf.append(
        #    np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        #)
        
        data_eulbias = np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{redshift:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
        
        data_lagbias = np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{redshift:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)

        b1l_idx = np.argmax(data_lagbias[:,1])
        b1e_idx = np.argmax(data_eulbias[:,1])

        b1l_kmax[ax_idx, o_idx] = data_lagbias[b1l_idx, 0]
        b1e_kmax[ax_idx, o_idx] = data_eulbias[b1e_idx, 0]        
        
        b1l_max[ax_idx, o_idx] = data_lagbias[b1l_idx, 1]
        b1e_max[ax_idx, o_idx] = data_eulbias[b1e_idx, 1]
        
        print("\t"+f"{b1l_max[ax_idx, o_idx]:.3f}"+"\t"+f"{b1e_max[ax_idx, o_idx]:.3f}")
    
        os.system('mv ./run.ini ./run_'+str(ax_idx)+'.ini')
        os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini')  

np.savetxt(rfpath+"plots/Figure_9_b1l_max.txt", b1l_max)
np.savetxt(rfpath+"plots/Figure_9_b1e_max.txt", b1e_max)
np.savetxt(rfpath+"plots/Figure_9_b1l_kmax.txt", b1l_kmax)
np.savetxt(rfpath+"plots/Figure_9_b1e_kmax.txt", b1e_kmax)

sns.set_style("white")

#Z = np.clip(np.nan_to_num(b1l_max), a_min=0., a_max=None)
Z = np.nan_to_num(b1l_max)
z_min = np.min(Z)
z_max = np.max(Z)

X, Y = np.meshgrid(m_ax, omega_ax/omega_cdm_LCDM)

for ax_idx, ax_val in enumerate(m_ax): 
    for o_idx, o_val in enumerate(omega_ax):
        if (Z[ax_idx, o_idx]==0.):
            Z[ax_idx, o_idx]=-1

#fig, ax = plt.subplots()
#c = ax.pcolormesh(X, Y, Z, cmap='RdBu', vmin=z_min, vmax=z_max)
#ax.set_title('pcolormesh')
## set the limits of the plot to the limits of the data
#ax.axis([x.min(), x.max(), y.min(), y.max()])
#fig.colorbar(c, ax=ax)
#plt.show()

#ax = sns.heatmap(Z)

fig, ax = plt.subplots(1,1, figsize=(15, 15))
hmap = ax.pcolormesh(X, Y, np.transpose(Z), vmin=-1, vmax=z_max)
ax.set_xscale('log')
ax.tick_params(axis='both', labelsize=25)
ax.set_xlabel("$M_\phi$ [eV]", fontsize=30)
ax.set_ylabel("$\Omega_\phi / \Omega_d$", fontsize=30)
cbar = plt.colorbar(hmap)
cbar.set_label(label="$\max(b^1_L)$", size=30)
cbar.ax.tick_params(labelsize=25) 

#fig, ax = plt.subplots(1, 1)
#ax.contourf(X, Y, np.transpose(Z), vmin=0., vmax=z_max)
#ax.set_xscale('log')

plt.savefig(rfpath+"plots/Figure_9.png")

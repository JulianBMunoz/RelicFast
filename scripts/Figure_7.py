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
    np.power(10., -25.), 
    np.power(10., -30.), 
    np.power(10., -32.),
])
omega_ax = 0.05*omega_cdm_LCDM
sum_massive_nu = 0.
redshifts = np.linspace(0., 10.0, 11) 
kmin = 1.0e-4
kmax = 1.0
Nk = 50

######################################################

omega_nu = sum_massive_nu/93.2
Mnu = sum_massive_nu # Units: eV

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = Mnu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = 0.70148
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

kfs = np.pi * np.sqrt(m_ax*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)

######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

data_tf = []
data_eulbias = []
data_lagbias = []
axion_background = []
axion_aosc = []

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
    for z_idx, z_val in enumerate(redshifts): 
        print("Running RelicFast + axionCAMB for m_ax = ", ax_val, ", z = ", z_val)
        
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
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{z_val:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                #"hubble = 0.701", "hubble = "+f'{h[ax_idx]:.6f}'
                "hubble = 0.701", "hubble = "+f'{h:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax[ax_idx]:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm[ax_idx]:.6e}'
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
        axion_background.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_background.dat"))
        axion_aosc.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_aosc.dat"))
        
        # Collect power spectrum for requested redshift 
        output_dirs = os.listdir(rfpath_boltzmannsuffix)
        output_dirs.remove('_params.ini')
        
        z_vals = np.array([])
        
        for str_idx, str_val in enumerate(output_dirs):
            if (str_val[0:15]=='_transfer_out_z'): 
                z_vals = np.append(
                    z_vals, 
                    np.float(str_val.split('_transfer_out_z')[1])
                )
                    
        z_vals = np.sort(z_vals)
        
        #PM = []
        #TF = []
        #for z_idx, z_val in enumerate(z_vals): 
        #    PM.append(np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(len(z_vals)-z_idx)+'.dat'))
        #    TF.append(np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_val:.3f}'))
        #
        #data_pm.append(PM)
        #data_tf.append(TF)
        
        data_eulbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_val:.2f}'
                +'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
        )
        data_lagbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_val:.2f}'
                +'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
        )
    
        os.system('mv ./run.ini ./run_'+str(ax_idx)+'.ini')
        os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini')  




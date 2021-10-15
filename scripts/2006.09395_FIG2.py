import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()


##################################
####  Setup Analysis         #####
##################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"

Mnu = 90.0e-3 # Units: eV
#Mnu = 0. # Units: eV
redshifts = np.linspace(0.65, 1.65, 11)[0:11]
kmin = 1.0e-4
kmax = 0.2
Nk = 50

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = Mnu/93.2

h = (0.701)*((omega_b_LCDM+omega_cdm_LCDM+omega_nu)/(omega_b_LCDM+omega_cdm_LCDM))

f_cdm = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

##################################
####  Run Boltzmann Solvers  #####
##################################

# Set solver to CAMB
camb_data_eulbias = []
camb_data_tf = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
        "define boltzmann_tag  _CLASS_", "define boltzmann_tag  _CAMB_"
    ).replace(
        "define boltzmann_tag  _AXIONCAMB_", "define boltzmann_tag  _CAMB_"
    )
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

# Calculate Bias using CAMB 
for z_idx, z_val in enumerate(redshifts): 
    print("Running RelicFast + CAMB for m_nu = ", Mnu, ", z = ", z_val)
    
    os.system('rm -r '+rfpath_outputsuffix)
    os.system('rm -r '+'Boltzmann_1')
    os.system('rm ./run.ini')
    
    # RUN RelicFast + CAMB for each neutrino mass choice 
    reading_file = open("2006.09395.ini", "r")
    new_file_content = ""
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(
            "mnu1 = 0.0", "mnu1 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{Mnu/3.:.2e}'
        ).replace(
            "tag_sterile_nu = 0", "tag_sterile_nu = 1"
        ).replace(
            "hubble=0.701", "hubble="+f'{h:.4f}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{z_val:.3f}'
        ).replace(
            "kbot= 1e-4", "kbot= "+f'{kmin:.3e}'
        ).replace(
            "ktop= 0.7", "ktop= "+f'{kmax:.3e}'
        ).replace(
            "N_klong = 1", "N_klong = "+f'{Nk:d}'
        )
        new_file_content += "    " + new_line +"\n"
    reading_file.close()
    writing_file = open("run.ini", "w")
    writing_file.write(new_file_content)
    writing_file.close()
    
    os.system('./relicfast run.ini')
    
    # Collect Eulerian bias for requested redshift 
    camb_data_eulbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_val:.2f}'+'_M13.00_Nk50.dat', skiprows=1)
    )
    camb_data_tf.append(
        np.loadtxt(rfpath_outputsuffix+'z'+f'{z_val:.2f}'+'TF_CAMB.dat', skiprows=1)
    )

    os.system('mv ./run.ini ./run_camb_'+str(z_idx)+'.ini')

# Set solver to CLASS
class_data_eulbias = []
class_data_tf = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
            "define boltzmann_tag  _CAMB_", "define boltzmann_tag  _CLASS_" 
    ).replace(
        "define boltzmann_tag  _AXIONCAMB_", "define boltzmann_tag  _CLASS_"
    )
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

# Calculate Bias using CLASS
for z_idx, z_val in enumerate(redshifts): 
    print("Running RelicFast + CLASS for m_nu = ", Mnu, ", z = ", z_val)
    
    os.system('rm -rf '+rfpath_outputsuffix)
    os.system('rm -rf '+'Boltzmann_0')
    os.system('rm ./run.ini')
    
    # RUN RelicFast + CLASS for each neutrino mass choice 
    reading_file = open("2006.09395.ini", "r")
    new_file_content = ""
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(
            "mnu1 = 0.0", "mnu1 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{Mnu/3.:.2e}'
        ).replace(
            "tag_sterile_nu = 0", "tag_sterile_nu = 1"
        ).replace(
            "hubble=0.701", "hubble="+f'{h:.4f}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{z_val:.3f}'
        ).replace(
            "kbot= 1e-4", "kbot= "+f'{kmin:.3e}'
        ).replace(
            "ktop= 0.7", "ktop= "+f'{kmax:.3e}'
        ).replace(
            "N_klong = 1", "N_klong = "+f'{Nk:d}'
        )
        new_file_content += "    " + new_line +"\n"
    reading_file.close()
    writing_file = open("run.ini", "w")
    writing_file.write(new_file_content)
    writing_file.close()
    
    os.system('./relicfast run.ini')
    
    # Collect Eulerian bias for requested redshift 
    class_data_eulbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_val:.2f}'+'_M13.00_Nk50.dat', skiprows=1)
    )
    class_data_tf.append(
        np.loadtxt(rfpath_outputsuffix+'z'+f'{z_val:.2f}'+'tk.dat', skiprows=1)
    )

    os.system('mv ./run.ini ./run_class_'+str(z_idx)+'.ini')

# Set solver to axionCAMB
axioncamb_data_eulbias = []
axioncamb_data_tf = []
axioncamb_data_cl = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
            "define boltzmann_tag  _CLASS_", "define boltzmann_tag  _AXIONCAMB_"
    ).replace(
        "define boltzmann_tag  _CAMB_", "define boltzmann_tag  _AXIONCAMB_"
    )
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

# Calculate Bias using axionCAMB 
for z_idx, z_val in enumerate(redshifts): 
    print("Running RelicFast + axionCAMB for m_nu = ", Mnu, ", z = ", z_val)
    
    os.system('rm -rf '+rfpath_outputsuffix)
    os.system('rm -rf '+'Boltzmann_2')
    os.system('rm ./run.ini')
    
    # RUN RelicFast + axionCAMB for each neutrino mass choice 
    reading_file = open("2006.09395.ini", "r")
    new_file_content = ""
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(
            "mnu1 = 0.0", "mnu1 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{Mnu/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{Mnu/3.:.2e}'
        ).replace(
            "tag_sterile_nu = 0", "tag_sterile_nu = 1"
        ).replace(
            "hubble=0.701", "hubble="+f'{h:.4f}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{z_val:.3f}'
        ).replace(
            "kbot= 1e-4", "kbot= "+f'{kmin:.3e}'
        ).replace(
            "ktop= 0.7", "ktop= "+f'{kmax:.3e}'
        ).replace(
            "N_klong = 1", "N_klong = "+f'{Nk:d}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = 1.0e-9"
        ).replace(
            "m_ax = 1.0e-22", "m_ax = 1.0e-22"
        )
        new_file_content += "    " + new_line +"\n"
    reading_file.close()
    writing_file = open("run.ini", "w")
    writing_file.write(new_file_content)
    writing_file.close()
    
    os.system('./relicfast run.ini')
    
    # Collect Eulerian bias for requested redshift 
    axioncamb_data_eulbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_val:.2f}'+'_M13.00_Nk50.dat', skiprows=1)
    )
    axioncamb_data_tf.append(
        np.loadtxt(rfpath_outputsuffix+'z'+f'{z_val:.2f}'+'TF_CAMB.dat', skiprows=1)
    )
    #axioncamb_data_cl.append(
    #    np.loadtxt('./Boltzmann_2/transfer_files_0/_scalCls.dat')
    #)
    

    os.system('mv ./run.ini ./run_axioncamb_'+str(z_idx)+'.ini')


##################################
####  Plot Results           #####
##################################

kplot = np.geomspace(10**-3.9, 0.1, 32)

plt.figure(figsize=(15, 7.5))
plt.xscale('log')

rgb = np.zeros(3)

for z_idx, z_val in enumerate(redshifts): 
    kvals = camb_data_eulbias[z_idx][:, 0]
    eulbiasvals = camb_data_eulbias[z_idx][:, 1]
    
    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
    eulbiasplot = eulbiasinterp(kplot)
    
    rgb[0] = 3./255.
    rgb[1] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    rgb[2] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    
    lines = [
        Line2D([0], [0], color='black', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='cyan', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='red', linewidth=3, linestyle='-')
    ]
    
    labels = [
        r'z = 1.65', 
        r'z = 0.65',
        r'$k_{\rm eq}$'
    ]
    
    plt.plot(kplot, eulbiasplot/eulbiasplot[0], color=tuple(rgb))

#plt.plot([0.7*0.015, 0.7*0.015], [0.985, 1.015], color='red')
plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=15)
plt.ylabel(r'$b_1(k)/b_1(k_{\rm ref})$', fontsize=15)
plt.title(r'CAMB, $ \Sigma m_\nu = $'+f'{Mnu:.2f}'+'meV', fontsize=15)
plt.legend(lines, labels, fontsize=15)
plt.grid(True, which='both', axis='both')
plt.saveplt('2006.09395_FIG2_CAMB.png')

kplot = np.geomspace(10**-3.9, 0.1, 32)

plt.figure(figsize=(15, 7.5))
plt.xscale('log')

rgb = np.zeros(3)

for z_idx, z_val in enumerate(redshifts): 
    kvals = class_data_eulbias[z_idx][:, 0]
    eulbiasvals = class_data_eulbias[z_idx][:, 1]
    
    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
    eulbiasplot = eulbiasinterp(kplot)
    
    rgb[0] = 3./255.
    rgb[1] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    rgb[2] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    
    lines = [
        Line2D([0], [0], color='black', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='cyan', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='red', linewidth=3, linestyle='-')
    ]
    
    labels = [
        r'z = 1.65', 
        r'z = 0.65',
        r'$k_{\rm eq}$'
    ]
    
    plt.plot(kplot, eulbiasplot/eulbiasplot[0], color=tuple(rgb))

#plt.plot([0.7*0.015, 0.7*0.015], [0.985, 1.015], color='red')
plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=15)
plt.ylabel(r'$b_1(k)/b_1(k_{\rm ref})$', fontsize=15)
plt.title(r'CLASS, $\Sigma m_\nu = $'+f'{Mnu:.1f}'+'meV', fontsize=15)
plt.legend(lines, labels, fontsize=15)
plt.grid(True, which='both', axis='both')
#plt.ylim((0.999, 1.010))
plt.saveplt('2006.09395_FIG2_CLASS.png')

kplot = np.geomspace(10**-3.9, 0.1, 32)

plt.figure(figsize=(15, 7.5))
plt.xscale('log')

rgb = np.zeros(3)

for z_idx, z_val in enumerate(redshifts): 
    kvals = axioncamb_data_eulbias[z_idx][:, 0]
    eulbiasvals = axioncamb_data_eulbias[z_idx][:, 1]
    
    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
    eulbiasplot = eulbiasinterp(kplot)
     
    rgb[0] = 3./255.
    rgb[1] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    rgb[2] = (len(redshifts)-1-z_idx) * (255/(len(redshifts)-1)) / 255.
    
    lines = [
        Line2D([0], [0], color='black', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='cyan', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='red', linewidth=3, linestyle='-')
    ]
    
    labels = [
        r'z = 1.65', 
        r'z = 0.65',
        r'$k_{\rm eq}$'
    ]
    
    plt.plot(kplot, eulbiasplot/eulbiasplot[0], color=tuple(rgb))
    
#plt.plot([0.7*0.015, 0.7*0.015], [0.985, 1.015], color='red')
plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=15)
plt.ylabel(r'$b_1(k)/b_1(k_{\rm ref})$', fontsize=15)
plt.title(r'axionCAMB, $ \Sigma m_\nu = $'+f'{Mnu:.2f}'+'meV', fontsize=15)
plt.legend(lines, labels, fontsize=15)
plt.grid(True, which='both', axis='both')
plt.saveplt('2006.09395_FIG2_AXIONCAMB.png')    
    
        
    
        
    
    

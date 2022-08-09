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
sns.set_style(style='white')

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
    10**-31,
    10**-32
])
omega_ax = np.array(len(m_ax)*[0.05*omega_cdm_LCDM])

sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-4
kmax = 0.5
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

h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

kfs = np.pi * np.sqrt(m_ax*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)

######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

data_pm = []
data_tf = []
data_eulbias = []
data_lagbias = []
axion_background = []

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
    print("Running RelicFast + axionCAMB for m_ax = ", ax_val)
    
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
            "hubble = 0.701", "hubble = "+f'{h[ax_idx]:.6f}'
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
        

    pm_idx = np.argmin(np.abs(z_vals - redshift))
    #pm_idx = np.argmin(np.abs(z_vals - 0.0))
    print('Requested/found redshift: ', redshift, z_vals[pm_idx])
    
    print('Loading: ')
    
    print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
    print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
    
    data_pm.append(
        np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
    )
    
    data_tf.append(
        np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
    )
    
    data_eulbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    )
    data_lagbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    )
    #print(data_lagbias[-1][0])

    os.system('mv ./run.ini ./run_'+str(ax_idx+8)+'.ini')
    os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx+8)+'.ini')  

kplot = np.geomspace(10**-3.9, 0.6, 100)


# eulerian bias plots
#plt.figure(figsize=(15, 10))
#plt.xscale('log')
#for ax_idx, ax_val in enumerate(m_ax): 
#    if ((ax_idx==0) or (ax_idx==(len(m_ax)-1))): 
#        continue 
#    
#    kvals = data_eulbias[ax_idx][:, 0]
#    eulbiasvals = data_eulbias[ax_idx][:, 1]
#    
#    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
#    eulbiasplot = eulbiasinterp(kplot)
#    
#    rgb = np.zeros(3)
#    rgb[0] = 3./255.
#    rgb[1] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
#    rgb[2] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
#    
#    yplot = eulbiasplot/eulbiasplot[0]
#    
#    plt.plot(kplot, yplot, label=r'$m_{\chi, i}=$'+f'{ax_val:.3e}'+r' eV', color=tuple(rgb))
#    plt.plot([kfs[ax_idx], kfs[ax_idx]], [min(yplot), max(yplot)], color=tuple(rgb), linestyle='dashed')
#    
#plt.plot([0.7*0.015, 0.7*0.015], [1., 1.01], color='red', label=r'$k_{eq}$')
#plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=30)
#plt.ylabel(r'$b_1(k)/b_1(k_{\rm ref})$', fontsize=30)
#plt.title(r'$\omega_\chi = 0.05\times\omega_{cdm},  ~z = 0.65, ~\Sigma M_\nu = 0$ eV', fontsize=30)
##plt.legend(lines, labels, fontsize=15)
#plt.legend(fontsize=25)
##plt.yscale('log')
#plt.grid(False, which='both', axis='both')

# lagrangian bias plots 
plt.figure(figsize=(15, 15))
plt.xscale('log')
for ax_idx, ax_val in enumerate(m_ax): 
    #if ((ax_idx==0) or (ax_idx==(len(m_ax)-1))):
    #    continue 
    
    kvals = data_lagbias[ax_idx][:, 0]
    lagbiasvals = data_lagbias[ax_idx][:, 1]
    
    lagbiasinterp = scipy.interpolate.interp1d(kvals, lagbiasvals)
    lagbiasplot = lagbiasinterp(kplot)
    
    rgb = np.zeros(3)
    rgb[0] = 3./255.
    rgb[1] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
    rgb[2] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
    
    yplot = lagbiasplot/lagbiasplot[0]
    
    plt.plot(
        kplot, 
        yplot, 
        label=r'$m_{\phi}= 10^{'+f'{np.log10(ax_val):.1f}'+r'}$ eV', 
        color=tuple(rgb), 
        linewidth=5.)
    plt.plot(
        [kfs[ax_idx], kfs[ax_idx]], 
        [min(yplot), max(yplot)], 
        color=tuple(rgb), 
        linestyle='dashed', 
        linewidth=5.)
    
plt.plot(
    [0.70148*0.015, 0.70148*0.015], 
    [1., 1.07], 
    color='red', 
    label=r'$k_{\rm eq}$',
    linewidth=5.)

plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=30)
plt.ylabel(r'$b_1^L(k)/b_1^L(k_{\rm ref})$', fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim((1.0e-4, 1.0e0))
#plt.title(r'$\omega_\chi = 0.05\times\omega_{cdm},  ~z = 0.65, ~\Sigma M_\nu = 0$ eV', fontsize=30)
#plt.legend(lines, labels, fontsize=15)
plt.legend(fontsize=25)
#plt.yscale('log')
plt.grid(False, which='both', axis='both')
plt.savefig(rfpath+"Figure_6.png") 


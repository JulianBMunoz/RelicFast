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

neutrino_masses = np.array([0.0e-3, 60.0e-3, 90.0e-3, 120.0e-3]) # Units: eV
redshifts = np.array([0.65, 1.15, 1.65])

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"

rfpath_outputsuffix = "output/result-0/"
outpath = "/Users/nicholasdeporzio/Desktop/"
reference_datapath = "/Users/nicholasdeporzio/Downloads/"

solver = "_CLASS_"

######################################################
####         INTERNAL 
######################################################

if solver=="_CLASS_":
    rfpath_boltzmannsuffix = "Boltzmann_0/transfer_files_0/"
elif solver=="_CAMB_":
    rfpath_boltzmannsuffix = "Boltzmann_1/transfer_files_0/"
elif solver=="_AXIONCAMB_":
    rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

data_pm = []
data_tf = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
            "#define boltzmann_tag  _CLASS_", "#define boltzmann_tag  "+solver
    ).replace(
        "#define boltzmann_tag  _CAMB_", "#define boltzmann_tag  "+solver
    ).replace(
        "#define boltzmann_tag  _AXIONCAMB_", "#define boltzmann_tag  "+solver
    )
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = neutrino_masses/93.2

h = np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu)/Omega_M_LCDM)

omega_cdm = omega_cdm_LCDM - omega_nu
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

os.chdir(rfpath)

for mnu_idx, mnu_val in enumerate(neutrino_masses): 
    print("Running RelicFast + ", solver, " for m_nu = ", mnu_val)
    
    os.system('rm -r '+rfpath_outputsuffix)
    os.system('rm -r '+'Boltzmann_0')
    os.system('rm -r '+'Boltzmann_1')
    os.system('rm -r '+'Boltzmann_2')
    os.system('rm ./run.ini')
    
    # RUN RelicFast + axionCAMB for each neutrino mass choice 
    reading_file = open("2006.09395.ini", "r")
    new_file_content = ""
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(
            "mnu1 = 0.0", "mnu1 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = 0.65"
        ).replace(
            "hubble = 0.701", "hubble = "+f'{h[mnu_idx]:.5f}'
        )
        
        if (mnu_val!=0.): 
            new_line = new_line.replace(
                "tag_sterile_nu = 0", "tag_sterile_nu = 1"
            )
            
        new_file_content += "    " + new_line +"\n"
    reading_file.close()
    writing_file = open("run.ini", "w")
    writing_file.write(new_file_content)
    writing_file.close()
    
    os.system('./relicfast run.ini')
    
    # Collect power spectrum for requested redshift 
    output_dirs = os.listdir(rfpath_boltzmannsuffix)
    if (solver!='_CLASS_'):
        output_dirs.remove('_params.ini')
    
    z_vals = np.array([])
    
     
    if (solver=='_CLASS_'):
        z_vals = np.loadtxt(rfpath+rfpath_boltzmannsuffix+'zlist.dat', skiprows=1)
    else: 
        for str_idx, str_val in enumerate(output_dirs):
            if (str_val[0:15]=='_transfer_out_z'): 
                z_vals = np.append(
                    z_vals, 
                    np.float(str_val.split('_transfer_out_z')[1])
                )
            
    if (solver!="_CLASS_"):      
        z_vals = np.sort(z_vals)
        
    for z_idx, z_val in enumerate(redshifts):
        if (solver=="_CLASS_"):
            pm_idx = z_vals[:,1][np.argmin(np.abs(z_vals[:,0] - z_val))]
            print('Requested/found redshift: ', z_val, z_vals[:,0][np.argmin(np.abs(z_vals[:,0] - z_val))])
        else: 
            pm_idx = np.argmin(np.abs(z_vals - z_val))
            print('Requested/found redshift: ', z_val, z_vals[pm_idx])
        
        print('Loading: ')
        
        if (solver=='_CLASS_'):
            print('\t '+rfpath_boltzmannsuffix+'z'+str(int(pm_idx+1))+'_tk.dat')
            data_tf.append(
                np.loadtxt(rfpath_boltzmannsuffix+'z'+str(int(pm_idx+1))+'_tk.dat')
            )
        else: 
            print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        
            data_pm.append(
                np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            )
        
            data_tf.append(
                np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
            )

    os.system('mv ./run.ini ./run_'+str(mnu_idx)+'.ini')
    

kplot = np.logspace(-4, 0, 41)

plt.figure(figsize=(15, 7.5))
plt.xscale('log')
for idx in range(len(data_tf)): 
    z_idx = idx%len(redshifts)
    z = redshifts[z_idx]
    m_idx = idx//len(redshifts)
    m = neutrino_masses[m_idx]
    
    kvals = data_tf[idx][:, 0]
    tcvals = data_tf[idx][:, 1]
    tbvals = data_tf[idx][:, 2]
    
    tcinterp = scipy.interpolate.interp1d(kvals, tcvals)
    tbinterp = scipy.interpolate.interp1d(kvals, tbvals)
    
    tcplot = tcinterp(kplot)
    tbplot = tbinterp(kplot)
    
    LCDM_kvals = data_tf[z_idx][:, 0]
    LCDM_tcvals = data_tf[z_idx][:, 1]
    LCDM_tbvals = data_tf[z_idx][:, 2]
    
    LCDM_tcinterp = scipy.interpolate.interp1d(LCDM_kvals, LCDM_tcvals)
    LCDM_tbinterp = scipy.interpolate.interp1d(LCDM_kvals, LCDM_tbvals)
    
    LCDM_tcplot = LCDM_tcinterp(kplot)
    LCDM_tbplot = LCDM_tbinterp(kplot)
    
    pcbtf = np.power(f_b[m_idx]*tbplot + f_cdm[m_idx]*tcplot, 2.)
    pcbtf_LCDM = np.power(f_b_LCDM*LCDM_tbplot + f_cdm_LCDM*LCDM_tcplot, 2.)
    
    tf_ratio_plot = (pcbtf - pcbtf_LCDM)/pcbtf_LCDM
    
    ref1 = np.loadtxt(reference_datapath+'2006.09395_Fig1_m_60mev_z_0.65.csv', delimiter=',')
    ref2 = np.loadtxt(reference_datapath+'2006.09395_Fig1_m_90mev_z_0.65.csv', delimiter=',')
    ref3 = np.loadtxt(reference_datapath+'2006.09395_Fig1_m_120mev_z_0.65.csv', delimiter=',')
    
    
    if m_idx==0: 
        col = 'black'
    elif m_idx==1: 
        col = 'red'
    elif m_idx==2: 
        col = 'cyan'
    elif m_idx==3: 
        col = 'blue'
    
    if z_idx==0: 
        ls='solid'
    elif z_idx==1: 
        ls='dashed'
    elif z_idx==2: 
        ls='dashdot'
        
    if idx>=len(redshifts):   
        plt.plot(kplot, tf_ratio_plot, color=col, linestyle=ls)
        
    plt.plot(ref1[:,0], ref1[:,1], color='black', linestyle='solid', linewidth=0.5)
    plt.plot(ref2[:,0], ref2[:,1], color='black', linestyle='solid', linewidth=0.5)
    plt.plot(ref3[:,0], ref3[:,1], color='black', linestyle='solid', linewidth=0.5)
    
    lines = [
        Line2D([0], [0], color='red', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='cyan', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='blue', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='grey', linewidth=3, linestyle='-'),
        Line2D([0], [0], color='grey', linewidth=3, linestyle='dashed'),
        Line2D([0], [0], color='grey', linewidth=3, linestyle='dashdot'),
        Line2D([0], [0], color='black', linewidth=3, linestyle='-')
    ]
    
    labels = [
        r'$\Sigma m_\nu = 60$meV', 
        r'$\Sigma m_\nu = 90$meV', 
        r'$\Sigma m_\nu = 120$meV',
        r'$z = 0.65$',
        r'$z = 1.00$',
        r'$z = 1.65$', 
        r'2009.09393, z=0.65'
    ]
    
    plt.xlabel(r'$k$', fontsize=15)
    plt.ylabel(r'$(P_{cb} - P_{cb, \Lambda} )/P_{cb, \Lambda}$', fontsize=15)
    plt.legend(lines, labels, fontsize=15)
    plt.title(solver, fontsize=15)
    plt.grid(True, which='both', axis='both')
    plt.savefig("2006.09395_FIG1.png") 
    
    

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

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
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
    np.power(10., -32),
####################################
    10**-25, 
    10**-26, 
    10**-27, 
    10**-28, 
    10**-29, 
    10**-30, 
    np.power(10., -31),
    np.power(10., -32),
])

shift = int(len(m_ax)/2)
tflookup = {
    'k/h' : 0,
    'cdm' : 1, 
    'baryon' : 2, 
    'photons' : 3, 
    'nu_massless' : 4, 
    'nu_massive' : 5, 
    'axion' : 6
}

omega_ax = np.append(
    np.array(int(len(m_ax)/2)*[0.05*omega_cdm_LCDM]),
    np.array(int(len(m_ax)/2)*[1.0e-9*omega_cdm_LCDM])
)
sum_massive_nu = 0.
redshift = 0.65
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

data_pm = []
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
    
    PM = []
    TF = []
    for z_idx, z_val in enumerate(z_vals): 
        PM.append(np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(z_idx+1)+'.dat'))
        TF.append(np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_val:.3f}'))
    
    data_pm.append(PM)
    data_tf.append(TF)
    
    #data_eulbias.append(
    #    np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    #)
    #data_lagbias.append(
    #    np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    #)
    #print(data_lagbias[-1][0])

    os.system('mv ./run.ini ./run_'+str(ax_idx)+'.ini')
    os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini')  

def Pprim(k): #k in units of Mpc^-1
    val = (1.
        *(2.*np.power(np.pi, 2.))
        *np.power(k, -3.)
        *(2.23e-9) #As
        *np.power(k/0.05, 0.9666-1.) #(k/kp)^(ns-1)
    )
    return val 

def Pij(k, Ti, Tj, fi, fj): 
    val = (
        (fi * fj)
        * 2.*(Ti * Tj)
        *Pprim(k)
    )
    return val 


Pij = np.vectorize(Pij)
Pprim = np.vectorize(Pprim)

# Compare CAMB/axionCAMB LCDM Transfers

ztxt='0.000'

plot_x = np.logspace(-4, 0, 101) #Units Mpc^-1

cambout = np.loadtxt("/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/CAMB/fortran/Boltzmann_2/transfer_files_0/_transfer_out_z"+ztxt)

# Transfers T(k) where T is unitless? and k is [Mpc^-1]
cambinterp_cdm = scipy.interpolate.interp1d(cambout[:,0]*h, cambout[:,1]*np.power(cambout[:,0]*h, 2.))
cambinterp_b = scipy.interpolate.interp1d(cambout[:,0]*h, cambout[:,2]*np.power(cambout[:,0]*h, 2.))
cambinterp_nu = scipy.interpolate.interp1d(cambout[:,0]*h, cambout[:,4]*np.power(cambout[:,0]*h, 2.))
cambinterp_total = scipy.interpolate.interp1d(cambout[:,0]*h, cambout[:,6]*np.power(cambout[:,0]*h, 2.))

pm_idx = np.argmin(np.abs(z_vals - float(ztxt)))
TF_lcdm = data_tf[0+shift][pm_idx]
TF_lcdm_k = TF_lcdm[:, tflookup['k/h']]*h # Careful of units 
TF_lcdm_cdm = TF_lcdm[:, tflookup['cdm']]*np.power(TF_lcdm_k, 2.)
TF_lcdm_b = TF_lcdm[:, tflookup['baryon']]*np.power(TF_lcdm_k, 2.)
TF_lcdm_x = TF_lcdm[:, tflookup['axion']]*np.power(TF_lcdm_k, 2.)
TF_lcdm_cdm_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_cdm)
TF_lcdm_b_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_b)
TF_lcdm_x_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_x)

plot_y1 = cambinterp_cdm(plot_x)
plot_y2 = TF_lcdm_cdm_interp(plot_x)
plt.figure(figsize=(15,7.5))
plt.plot(plot_x, plot_y1, label='CAMB')
plt.plot(plot_x, plot_y2, label='axionCAMB', linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=16)
plt.xlabel("k", fontsize=16)
plt.ylabel("$T_{cdm}$", fontsize=16)

plot_y1 = cambinterp_b(plot_x)
plot_y2 = TF_lcdm_b_interp(plot_x)
plt.figure(figsize=(15,7.5))
plt.plot(plot_x, plot_y1, label='CAMB')
plt.plot(plot_x, plot_y2, label='axionCAMB', linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=16)
plt.xlabel("k", fontsize=16)
plt.ylabel("$T_{b}$", fontsize=16)

fcdm = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
fb = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
CAMB_Pmm_lcdm_reco_bcdm = (
    Pprim(plot_x)
    * (
        np.power(fb, 2.)*np.power(cambinterp_b(plot_x), 2.)
        + 2.*(fb*fcdm)*(cambinterp_b(plot_x)*cambinterp_cdm(plot_x))
        + np.power(fcdm, 2.)*np.power(cambinterp_cdm(plot_x), 2.)
    )
)
CAMB_Pmm_lcdm_reco_tot = (
    Pprim(plot_x)
    * np.power(cambinterp_total(plot_x), 2.)
)

cambout_pmm = np.loadtxt("/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/CAMB/fortran/Boltzmann_2/transfer_files_0/_matterpower_100.dat")
cambinterp_pmm = scipy.interpolate.interp1d(cambout_pmm[:,0]*h, cambout_pmm[:,1]*np.power(h, -3.))
plot_y1 = cambinterp_pmm(plot_x)
plot_y2 = CAMB_Pmm_lcdm_reco_bcdm
plot_y3 = CAMB_Pmm_lcdm_reco_tot
plt.figure(figsize=(15,7.5))
plt.plot(plot_x, plot_y1, label='CAMB Pmm Direct')
plt.plot(plot_x, plot_y2, label='CAMB Pmm Reco b+CDM', linestyle='dashed')
plt.plot(plot_x, plot_y3, label='CAMB Pmm Reco Total', linestyle='dotted')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=16)
plt.xlabel("k", fontsize=16)
plt.ylabel("$P_{mm}$", fontsize=16)
plt.savefig(rfpath+"plots/match_camb_axioncamb_pmm_1.png") 

plt.figure(figsize=(15,7.5))
plt.plot(plot_x, plot_y1/plot_y3, label='CAMB Pmm Direct/CAMB Pmm Reco')
plt.xscale('log')
#plt.yscale('log')
plt.legend(fontsize=16)
plt.xlabel("k", fontsize=16)
plt.ylabel("$P_{mm, direct}/P_{mm, reco}$", fontsize=16)
plt.savefig(rfpath+"plots/match_camb_axioncamb_pmm_2.png")

#z_table = np.array(z_vals)
#k_table = np.geomspace(1.0e-4, 2., 100)
#
#k_ref = 1.0
#
#z_plot = np.linspace(min(z_vals), max(z_vals), 100)
#colors=sns.color_palette("flare", int(len(m_ax)/2))
#
#plt.figure(figsize=(15,7.5))
#for ax_idx, ax_val in enumerate(m_ax[0:8]):   
#    tf_vals = np.zeros(len(z_table))
#    tf_lcdm_vals = np.zeros(len(z_table))
#    
#    for z_idx, z_val in enumerate(z_table): 
#        TF = data_tf[ax_idx][z_idx]
#        TF_lcdm = data_tf[ax_idx+8][z_idx]
#        
#        tf_vals[z_idx] = scipy.interpolate.interp1d(TF[:, tflookup["k/h"]], TF[:, tflookup[species]])(k_ref)
#        tf_lcdm_vals[z_idx] = scipy.interpolate.interp1d(TF_lcdm[:, tflookup["k/h"]], TF_lcdm[:, tflookup[species]])(k_ref)
#        
#    tf_plot = (tf_vals-tf_lcdm_vals)/tf_lcdm_vals
#    
#    plt.plot(z_table, tf_plot, label="$m_\chi = $"+f'{ax_val:.2e}', color=colors[ax_idx])
#
#plt.xlabel('z')
#plt.ylabel('$\Delta T/T_{m\Lambda CDM}$')
#plt.legend()
#plt.title('T_'+species+', $\omega_\chi = 0.05 \omega_{cdm}$, $k = $'+f'{k_ref:.2e}')
#plt.show()


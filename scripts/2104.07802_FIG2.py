import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()

relicfast_data_eulbias = []
relicfast_data_lagbias = []
axioncamb_data_tf = []
axioncamb_data_pk = []
relicfast_data_tf = []
relicfast_data_pk = []

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
rfpath_outputsuffix = "output/result-0/"
rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"
outpath = "/Users/nicholasdeporzio/Downloads/"

M_nu = 0. # Units: eV
redshift = 0.7
kmin = 5.0e-5
kmax = 1.5
Nk = 50
M_halo = 1.0e13

h_lcdm = 0.67
omega_nu = M_nu/93.2
omega_cdm_LCDM = 0.12
omega_b_LCDM = 0.022
Omega_M_LCDM = (omega_b_LCDM + omega_cdm_LCDM)/np.power(h_lcdm, 2.)

Omega_ax = np.array([1e-9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08])*omega_cdm_LCDM/np.power(h_lcdm, 2.)
omega_ax = Omega_ax*np.power(h_lcdm, 2.)
omega_cdm = omega_cdm_LCDM - omega_ax


M_ax = [1.0e-30, 1.0e-28, 1.0e-26]

# Set solver to axionCAMB
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
for m_idx, m_val in enumerate(M_ax):
    for o_idx, o_val in enumerate(omega_ax):
        print("Running RelicFast + axionCAMB for m_ax = ", M_ax[m_idx])
        
        os.system('rm -rf '+rfpath_outputsuffix)
        os.system('rm -rf '+'Boltzmann_2')
        os.system('rm ./run.ini')
        
        # RUN RelicFast + axionCAMB for each neutrino mass choice 
        reading_file = open("1805.11623.ini", "r")
        new_file_content = ""
        for line in reading_file:
            stripped_line = line.strip()
            new_line = stripped_line.replace(
                "mnu1 = 0.0", "mnu1 = "+f'{M_nu:.2e}'
            ).replace(
                "mnu2 = 0.0", "mnu2 = "+f'{M_nu:.2e}'
            ).replace(
                "m_SN = 0.02", "m_SN = "+f'{M_nu:.2e}'
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h_lcdm:.4f}'
            ).replace(
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "kbot = 1e-4", "kbot = "+f'{kmin:.3e}'
            ).replace(
                "ktop = 0.7", "ktop = "+f'{kmax:.3e}'
            ).replace(
                "N_klong = 1", "N_klong = "+f'{Nk:d}'
            ).replace(
                "omegab = 0.02226", "omegab = "+f'{omega_b_LCDM:.3f}'
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm[o_idx]:.3f}'
            ).replace(
                "Mhalo_min = 1e13", "Mhalo_min = "+f'{M_halo/h_lcdm:.3e}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{o_val:.3e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{m_val:.3e}'
            )
            if (M_nu!=0.0): 
                new_line = new_line.replace(
                    "tag_sterile_nu = 0", "tag_sterile_nu = 1"
                )
            new_file_content += "    " + new_line +"\n"
        reading_file.close()
        writing_file = open("run.ini", "w")
        writing_file.write(new_file_content)
        writing_file.close()
        
        os.system('./relicfast run.ini')
        
        relicfast_data_eulbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{redshift:.2f}'+'_M'
                       +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                       +'_Nk50.dat', skiprows=1)
        )
        relicfast_data_lagbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{redshift:.2f}'+'_M'
                       +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                       +'_Nk50.dat', skiprows=1)
        )
        relicfast_data_tf.append(
            np.loadtxt(rfpath_outputsuffix+'z'+f'{redshift:.2f}'+'TF_CAMB.dat', skiprows=1)
        )
        relicfast_data_pk.append(
            np.loadtxt(rfpath_outputsuffix+'power_spectra_'+'z'+f'{redshift:.2f}'+'_M'
                       +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                       +'_Nk50.dat', skiprows=1)
        )

        #############
        # Collect axioncamb power spectrum for requested redshift 
        #output_dirs = os.listdir(rfpath_boltzmannsuffix)
        #output_dirs.remove('_params.ini')
        #
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
        #    
    
        #pm_idx = np.argmin(np.abs(z_vals - redshift))
        #print('Requested/found redshift: ', redshift, z_vals[pm_idx])
        #
        #print('Loading: ')
        #
        #print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        #print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        #
        #axioncamb_data_pk.append(
        #    np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        #)
        #
        #axioncamb_data_tf.append(
        #    np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        #)

        #############
        
        os.system('mv ./run.ini ./run_axioncamb_'+str(m_idx)+'.ini')


kplot = np.geomspace(10**-3, 10**0.3, 34) #Units: h Mpc^-1

colors = sns.color_palette("magma", len(omega_ax))

for m_idx, m_val in enumerate(M_ax):
    plt.figure(figsize=(15, 7.5))
    ref = np.loadtxt(rfpath+"scripts/2104.07802_FIG2_REF_"+str(m_idx)+".csv", delimiter=",")
    
    for o_idx, o_val in enumerate(omega_ax): 
        idx = len(omega_ax)*m_idx + o_idx
        
        if o_idx==0: 
            k_lcdm_vals = relicfast_data_pk[idx][:, 0]/h_lcdm
            pm_lcdm_vals = relicfast_data_pk[idx][:, 1]
            pm_lcdm_interp = scipy.interpolate.interp1d(k_lcdm_vals, pm_lcdm_vals)
            pm_lcdm_plot = pm_lcdm_interp(kplot)
        
        kvals = relicfast_data_pk[idx][:, 0]/h_lcdm # h/Mpc
        pmvals = relicfast_data_pk[idx][:, 1] # Pmm
    
        pminterp = scipy.interpolate.interp1d(kvals, pmvals)
        pmplot = pminterp(kplot)
        
        plt.plot(
            kplot, 
            pmplot/pm_lcdm_plot, 
            label=r'RelicFast+axionCAMB, $\Omega_{ax}/\Omega_{d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}', 
            color=colors[o_idx]
        )
        
        if (o_idx==(len(omega_ax)-1)):
            plt.scatter(ref[:,0][::10], ref[:,o_idx+1][::10], marker='.', color="black", label="2104.07802")
        else: 
            plt.scatter(ref[:,0][::10], ref[:,o_idx+1][::10], marker='.', color="black")
    
    plt.xscale('log')
    plt.xlabel(r'$k ~[h ~{\rm Mpc}^{-1}]$', fontsize=15)
    plt.ylabel(r'$P_m/P_{m, ~LCDM}$', fontsize=15)
    plt.title(r'2104.07802 FIG 2, axionCAMB, $M_\chi = $'+f'{m_val:.3e}', fontsize=15)
    plt.legend(fontsize=15)
    plt.grid(True, which='both', axis='both')
    plt.savefig("2104.07802_FIG2"+str(m_idx)+".png")    

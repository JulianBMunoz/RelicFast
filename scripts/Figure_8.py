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

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"
outpath = "/Users/nicholasdeporzio/Desktop/"

sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-5
kmax = 0.5
Nk = 50

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

######################################################

omega_nu = sum_massive_nu/93.2

kref = np.power(10., -4.) 
omega_ax = np.array([omega_cdm_LCDM*1.0e-12, omega_cdm_LCDM*0.05])
m_ax = np.array([
    np.power(10., -32.),
    np.power(10., -28.),
    np.power(10., -26.)
])

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = sum_massive_nu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = 0.70148
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

data_pm = []
data_tf = []
data_eulbias = []
data_lagbias = []
axion_background = []
relicfast_pss = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    if "#define length_transfer_axioncamb " in stripped_line:
        expected_axioncamb_output_lines = int(''.join(filter(str.isdigit, stripped_line)))
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
for m_idx, m_val in enumerate(m_ax): 
    for oax_idx, oax_val in enumerate(omega_ax): 
        print("Running RelicFast + axionCAMB for omega_ax = ", oax_val)
        
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
                "mnu1 = 0.0", "mnu1 = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "mnu2 = 0.0", "mnu2 = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "m_SN = 0.02", "m_SN = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{oax_val:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm[oax_idx]:.6e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{m_val:.6e}'
            ).replace(
                "kbot = 1e-4", "kbot = "+f'{kmin:.3e}'
            )
            
            
            if (sum_massive_nu!=0.): 
                new_line = new_line.replace(
                    "tag_sterile_nu = 0", "tag_sterile_nu = 1"
                )
                
            new_file_content += "    " + new_line +"\n"
        reading_file.close()
        writing_file = open("run.ini", "w")
        writing_file.write(new_file_content)
        writing_file.close()
        
        lines_match_flag=False
        while lines_match_flag==False:
            os.system('./relicfast run.ini')
            with open(rfpath+"/Boltzmann_2/transfer_files_0/_transfer_out_z200.000", 'r') as fp:
                ac_lines_out = sum(1 for line in fp)
            fp.close()
    
            if ac_lines_out==expected_axioncamb_output_lines:
                lines_match_flag=True
            else:
                print("axionCAMB output doesn't match RelicFast compile parameters. Recompiling: "
                    +str(expected_axioncamb_output_lines)+"-->"+str(ac_lines_out))
                os.chdir(rfpath+'/include/')
                reading_file = open("common.h", "r")
                new_file_content = ""
                for line in reading_file:
                    stripped_line = line.strip()
                    new_line = stripped_line.replace(
                            "#define length_transfer_axioncamb "+str(expected_axioncamb_output_lines),
                            "#define length_transfer_axioncamb "+str(ac_lines_out)
                    )
                    new_file_content += "    " + new_line +"\n"
                reading_file.close()
                writing_file = open("common.h", "w")
                writing_file.write(new_file_content)
                writing_file.close()
                os.chdir(rfpath)
                os.system('make')
    
                expected_axioncamb_output_lines = int(ac_lines_out)
        
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
    
        relicfast_pss.append(np.loadtxt(
            rfpath_outputsuffix+'power_spectra_z'+f"{redshift:.2f}"+"_M13.00_Nk"+f"{Nk:d}"+".dat"
            , skiprows=1)
        )
    
        os.system('mv ./run.ini ./run_'+str(oax_idx)+'.ini')
        os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(oax_idx)+'.ini')    

#####################################

colors = sns.color_palette('magma', len(m_ax))
kplot = np.geomspace(10**-4.0, 0.1, 100)

fig, ax = plt.subplots(1, 1, figsize=(15., 15.))
for m_idx, m_val in enumerate(m_ax): 
    for oax_idx, oax_val in enumerate(omega_ax): 
        k_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 0]
        pmm_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 1]
        phh_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 3]
        
        pmm_interp = scipy.interpolate.interp1d(k_vals, pmm_vals)
        phh_interp = scipy.interpolate.interp1d(k_vals, phh_vals)
    
        if oax_idx==0:
            pmm_LCDM = pmm_interp(kplot)
            phh_LCDM = phh_interp(kplot)
            pmm_LCDM_ref = pmm_interp(kref)
            phh_LCDM_ref = phh_interp(kref) 
            print(pmm_LCDM)
            print(phh_LCDM)
    
        else:
            pmm = pmm_interp(kplot)
            phh = phh_interp(kplot) 
            Rmm = pmm/pmm_LCDM
            Rhh = phh/phh_LCDM
            Rmm_ref = pmm_interp(kref)/pmm_LCDM_ref
            Rhh_ref = phh_interp(kref)/phh_LCDM_ref
    
            yplot1 = Rmm/Rmm_ref
            yplot2 = Rhh/Rhh_ref    
     
            ax.plot(#Matter
                kplot,
                yplot1, 
                color=colors[m_idx], 
                linewidth=5., 
                linestyle='dashed'
            )
            ax.plot(#Halo 
                kplot,
                yplot2, 
                label=(r"$m_\phi = 10^{"+f"{np.log10(m_val):.0f}"+r"}$ eV"), 
                color=colors[m_idx], 
                linewidth=5.,
                linestyle='solid'
            )
#ax.plot([min(kplot), max(kplot)], [1., 1.], color='black', linewidth=5.)
ax.set_xlim((min(kplot), max(kplot)))
#ax.plot([0.7*0.015, 0.7*0.015], [0.999, 1.008], color='red', label=r'$k_{eq}$')
#ax.plot([0.024, 0.024], [0.999, 1.008], color='blue', label=r'$k_{*}$')
ax.set_xscale('log')
ax.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=40)
ax.set_ylabel(r'$R(k)/R(k_{\rm ref})$', fontsize=40)
ax.set_ylim((0.7, 1.002))
ax.tick_params(axis='both', labelsize=30)
ax.legend(fontsize=30)
ax.grid(False)
plt.savefig(rfpath+"plots/Figure_8.png")

#fig, ax = plt.subplots(1, 1, figsize=(20, 20))
#for oax_idx, oax_val in enumerate(omega_ax): 
#    k_vals = relicfast_pss[oax_idx][:, 0]
#    pmm_vals = relicfast_pss[oax_idx][:, 1]
#    phh_vals = relicfast_pss[oax_idx][:, 3]
#    
#    pmm_interp = scipy.interpolate.interp1d(k_vals, pmm_vals)
#    phh_interp = scipy.interpolate.interp1d(k_vals, phh_vals)
#
#    if oax_idx==0:
#        pmm_LCDM = pmm_interp(kplot)
#        phh_LCDM = phh_interp(kplot) 
#        print(pmm_LCDM)
#        print(phh_LCDM)
#
#    else:
#        pmm = pmm_interp(kplot)
#        phh = phh_interp(kplot) 
#        print(pmm)
#        print(phh)
#        yplot1 = pmm/pmm_LCDM
#        yplot2 = phh/phh_LCDM    
#
#        ax.plot(
#            kplot,
#            yplot2, 
#            label=r"$R_h(k)$", 
#            color=colors[oax_idx-1], 
#            linewidth=5.,
#            linestyle='dashed'
#        )
##ax.plot([min(kplot), max(kplot)], [1., 1.], color='black', linewidth=5.)
#ax.set_xlim((min(kplot), max(kplot)))
##ax.plot([0.7*0.015, 0.7*0.015], [0.999, 1.008], color='red', label=r'$k_{eq}$')
##ax.plot([0.024, 0.024], [0.999, 1.008], color='blue', label=r'$k_{*}$')
#ax.set_xscale('log')
#ax.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$', fontsize=40)
#ax.set_ylabel(r'$R(k)$', fontsize=40)
#ax.tick_params(axis='both', labelsize=30)
#ax.legend(fontsize=30)
#ax.grid(False)
#ax.set_title(r"$m_\phi = 10^{"+f"{np.log10(m_ax):.1f}"+r"}$ eV", fontsize=40)
#plt.savefig(rfpath+"plots/Figure_8b_logmax"+f"{np.log10(m_ax):.1f}"+".png")

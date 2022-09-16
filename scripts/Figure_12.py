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
from matplotlib import rc
import matplotlib.font_manager
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
matplotlib.rcParams.update({
    "font.weight" : "bold",
    "font.size" : 60,
    "axes.labelsize" : 110,
    "axes.labelpad" : 8.0,  
    "xtick.labelsize" : 60, 
    "ytick.labelsize" : 60, 
    "legend.fontsize" : 60, 
    "figure.dpi" : 300, 
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True, 
    "savefig.pad_inches" : 0.1

})
######################################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"

sum_massive_nu = 0.060
redshift = 0.65
kmin = 1.0e-4
kmax = 0.5
Nk = 50

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

######################################################

omega_nu = sum_massive_nu/93.2
Mnu = np.array([0.0, sum_massive_nu]) # Units: eV

#omega_ax = np.array([3.5*omega_nu, omega_cdm_LCDM*1.0e-12])
omega_ax = np.array([omega_nu, omega_cdm_LCDM*1.0e-12])

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = Mnu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = 0.70148
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

m_ax = ( #double check this
    np.power(0.08, 2.)
    * np.power(np.pi, -2.)
    * np.power(1+redshift, -2)
    * np.power(1000.*sum_massive_nu/3./100., 2.)
    * np.power(h, 2.)
    * np.power(10., 33.)
    * np.power(Omega_M_LCDM, -0.5)
    * np.power(1.+redshift, -1.5)
    * np.power(1.56, -2.)
    * np.power(10., -58)
)
#m_ax = 1.0e-30
print((
    "For single neutrino mass: "
    +f'{sum_massive_nu/3.:.3e}'
    +"\n, axion mass with same suppresion scale: "
    +f'{m_ax:.3e}'
    +"\n, and abundance little omega_ax = "
    f'{omega_ax[0]:.3e}'
))

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
for mnu_idx, mnu_val in enumerate(Mnu): 
    print("Running RelicFast + axionCAMB for m_nu = ", mnu_val)
    
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
            "mnu1 = 0.0", "mnu1 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
        ).replace(
            "z_collapse_top = 1.4", "z_collapse_top = 1.65"
        ).replace(
            "N_zcoll = 1", "N_zcoll = 1"
        ).replace(
            "hubble = 0.701", "hubble = "+f'{h:.6f}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax[mnu_idx]:.6e}'
        ).replace(
            "N_klong = 1", "N_klong = "+str(Nk)
        ).replace(
            "omegac = 0.11271", "omegac = "+f'{omega_cdm[mnu_idx]:.6e}'
        ).replace(
            "m_ax = 1.0e-22", "m_ax = "+f'{m_ax:.6e}'
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

    os.system('mv ./run.ini ./run_'+str(mnu_idx)+'.ini')
    os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(mnu_idx)+'.ini')    





#####################################

kplot = np.geomspace(10**-3.9, 0.1, 32)

fig, ax1 = plt.subplots(1, 1)
#fig.subplots_adjust(hspace=0)
for mnu_idx, mnu_val in enumerate(Mnu): 
    kvals = data_eulbias[mnu_idx][:, 0]
    eulbiasvals = data_eulbias[mnu_idx][:, 1]
    lagbiasvals = data_lagbias[mnu_idx][:, 1]
    
    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
    lagbiasinterp = scipy.interpolate.interp1d(kvals, lagbiasvals)

    eulbiasplot = eulbiasinterp(kplot)
    lagbiasplot = lagbiasinterp(kplot)
    
    rgb = np.zeros(3)
    rgb[0] = 3./255.
    rgb[1] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.
    rgb[2] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.

    if mnu_idx==0:
        labeltxt = r"$\phi ~+~ \rm{Massless Neutrinos}$"
    elif mnu_idx==1:
        labeltxt = r"${\rm Massive Neutrinos}$"
 
    ax1.plot(
        kplot, 
        lagbiasplot/lagbiasplot[0], 
        label=labeltxt, 
        color=tuple(rgb), 
        linewidth=5.
    )
ax1.set_xlim((min(kplot), max(kplot)))
#ax.plot([0.7*0.015, 0.7*0.015], [0.999, 1.008], color='red', label=r'$k_{eq}$')
#ax.plot([0.024, 0.024], [0.999, 1.008], color='blue', label=r'$k_{*}$')
ax1.set_xscale('log')
ax1.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$')
ax1.set_ylabel(r'$b_1^L(k)/b_1^L(k_{\rm ref})$')
ax1.tick_params(axis='both')
ax1.legend()
ax1.grid(False)
#plt.savefig(rfpath+"plots/Figure_10.png")

#for mnu_idx, mnu_val in enumerate(Mnu): 
#    kvals = data_eulbias[mnu_idx][:, 0]
#    eulbiasvals = data_eulbias[mnu_idx][:, 1]
#    lagbiasvals = data_lagbias[mnu_idx][:, 1]
#    
#    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
#    lagbiasinterp = scipy.interpolate.interp1d(kvals, lagbiasvals)
#
#    eulbiasplot = eulbiasinterp(kplot)
#    lagbiasplot = lagbiasinterp(kplot)
#    
#    rgb = np.zeros(3)
#    rgb[0] = 3./255.
#    rgb[1] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.
#    rgb[2] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.
#
#    if mnu_idx==0:
#        labeltxt = r"$\phi$ + Massless Neutrinos"
#    elif mnu_idx==1:
#        labeltxt = r"Massive Neutrinos"
# 
#    ax2.plot(
#        kplot, 
#        lagbiasplot, 
#        label=r'$m_{\nu, i}=$'+f'{1000.*mnu_val/3.:.0f}'+r' meV', 
#        color=tuple(rgb), 
#        linewidth=5.
#    )





sum_massive_nu = 1.0
redshift = 0.65
kmin = 1.0e-4
kmax = 0.5
Nk = 50

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

######################################################

omega_nu = sum_massive_nu/93.2
Mnu = np.array([0.0, sum_massive_nu]) # Units: eV

#omega_ax = np.array([3.5*omega_nu, omega_cdm_LCDM*1.0e-12])
omega_ax = np.array([omega_nu, omega_cdm_LCDM*1.0e-12])

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = Mnu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = 0.70148
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

m_ax = ( #double check this
    np.power(0.08, 2.)
    * np.power(np.pi, -2.)
    * np.power(1+redshift, -2)
    * np.power(1000.*sum_massive_nu/3./100., 2.)
    * np.power(h, 2.)
    * np.power(10., 33.)
    * np.power(Omega_M_LCDM, -0.5)
    * np.power(1.+redshift, -1.5)
    * np.power(1.56, -2.)
    * np.power(10., -58)
)
#m_ax = 1.0e-30
print((
    "For single neutrino mass: "
    +f'{sum_massive_nu/3.:.3e}'
    +"\n, axion mass with same suppresion scale: "
    +f'{m_ax:.3e}'
    +"\n, and abundance little omega_ax = "
    f'{omega_ax[0]:.3e}'
))

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
for mnu_idx, mnu_val in enumerate(Mnu): 
    print("Running RelicFast + axionCAMB for m_nu = ", mnu_val)
    
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
            "mnu1 = 0.0", "mnu1 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{mnu_val/3.:.2e}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
        ).replace(
            "z_collapse_top = 1.4", "z_collapse_top = 1.65"
        ).replace(
            "N_zcoll = 1", "N_zcoll = 1"
        ).replace(
            "hubble = 0.701", "hubble = "+f'{h:.6f}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax[mnu_idx]:.6e}'
        ).replace(
            "N_klong = 1", "N_klong = "+str(Nk)
        ).replace(
            "omegac = 0.11271", "omegac = "+f'{omega_cdm[mnu_idx]:.6e}'
        ).replace(
            "m_ax = 1.0e-22", "m_ax = "+f'{m_ax:.6e}'
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

    os.system('mv ./run.ini ./run_'+str(mnu_idx)+'.ini')
    os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(mnu_idx)+'.ini')    





#####################################

for mnu_idx, mnu_val in enumerate(Mnu): 
    kvals = data_eulbias[mnu_idx][:, 0]
    eulbiasvals = data_eulbias[mnu_idx][:, 1]
    lagbiasvals = data_lagbias[mnu_idx][:, 1]
    
    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
    lagbiasinterp = scipy.interpolate.interp1d(kvals, lagbiasvals)

    eulbiasplot = eulbiasinterp(kplot)
    lagbiasplot = lagbiasinterp(kplot)
    
    rgb = np.zeros(3)
    rgb[0] = 3./255.
    rgb[1] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.
    rgb[2] = (len(Mnu)-1-mnu_idx) * (255/(len(Mnu)-1)) / 255.

    if mnu_idx==0:
        labeltxt = r"$\phi {\rm ~+~ Massless Neutrinos}$"
    elif mnu_idx==1:
        labeltxt = r"${\rm Massive Neutrinos}$"
 
    ax1.plot(
        kplot, 
        lagbiasplot/lagbiasplot[0], 
        #label=labeltxt, 
        linestyle='dashed',
        color=tuple(rgb), 
        linewidth=5.
    )

plt.savefig(rfpath+"plots/Figure_12.png")

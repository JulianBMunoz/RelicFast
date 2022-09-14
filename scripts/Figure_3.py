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
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

axion_rho_of_a = []

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
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

Omega_ax = 0.03*omega_cdm_LCDM/np.power(h_lcdm, 2.)
omega_ax = Omega_ax*np.power(h_lcdm, 2.)

omega_cdm = omega_cdm_LCDM - omega_ax

M_ax = np.logspace(-32, -22, 3)

# Set solver to axionCAMB
os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    if "#define length_transfer_axioncamb " in stripped_line:
        expected_axioncamb_output_lines = int(''.join(filter(str.isdigit, stripped_line)))
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

# Calculate axion redshift at fixed mass, for various abundance w/ axionCAMB
for m_idx, m_val in enumerate(M_ax):
    print("Running RelicFast + axionCAMB for m_ax = "+f'{m_val:.2e}')
    
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
            "omegac = 0.11271", "omegac = "+f'{omega_cdm:.3f}'
        ).replace(
            "Mhalo_min = 1e13", "Mhalo_min = "+f'{M_halo/h_lcdm:.3e}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax:.3e}'
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

    axion_rho_of_a.append(
        np.loadtxt(outpath+'axion_grhoax_internal.dat')
    )

        
    os.system('mv ./run.ini ./run_fixed_mass_'+str(m_idx)+'.ini')

textvars = []
colors = sns.color_palette("magma", len(M_ax))
plt.figure(figsize=(20, 15))
for m_idx, m_val in enumerate(M_ax): 
    avals = axion_rho_of_a[m_idx][:, 0] 
    rhovals = axion_rho_of_a[m_idx][:, 1] 

    rhointerp = scipy.interpolate.interp1d(avals, rhovals)
    
    a_osc = np.loadtxt(outpath+"axion_aosc.dat")
    rho_osc = rhointerp(a_osc)
    a_dm = np.geomspace(avals[-1], 10**0, 50)
    #rho_dm = rho_osc*np.power(a_osc/a_dm, 3)
    rho_dm = rhovals[-1]*np.power(avals[-1]/a_dm, 3)

    plt.plot(
        np.linspace(np.log10(avals[0]), np.log10(avals[-1])+0.5, 10), 
        np.log10(10*[rhovals[0]]), 
        color=colors[m_idx],
        linewidth=5.0, 
        linestyle='dashed'
    )   
    plt.plot(
        np.log10(a_dm), 
        np.log10(rho_dm), 
        color=colors[m_idx],
        linewidth=5.0, 
        linestyle='dotted'
    )

    plt.plot(
        np.log10(avals), 
        np.log10(rhovals), 
        label=(r'$M_\phi = 10^{'+f'{int(np.log10(m_val)):d}'+r'}$ eV'), 
        linewidth=5.0, 
        color=colors[m_idx]
    )


plt.xscale('linear')
plt.yscale('linear')
plt.xlabel(r'$\log({a})$', fontsize=40)
plt.ylabel(r'$\log(\omega_\phi)$', fontsize=40)
plt.xlim((-7, 0.1))
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
#plt.title(r'$M_\chi = $'+f'{M_ax_fixed:.3e}', fontsize=15)
plt.legend(fontsize=30, loc='upper right')
plt.grid(False, which='both', axis='both')
plt.savefig(rfpath+"plots/Figure_3.png") 


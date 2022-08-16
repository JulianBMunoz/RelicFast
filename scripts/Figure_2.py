import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()
sns.set_style(style='white')

axion_rho_of_a = []

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

Omega_ax = np.array([0.01, 0.03, 0.05])*omega_cdm_LCDM/np.power(h_lcdm, 2.)
omega_ax = Omega_ax*np.power(h_lcdm, 2.)

omega_cdm = omega_cdm_LCDM - omega_ax

M_ax_fixed = np.power(10., -26)

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

# Calculate axion redshift at fixed mass, for various abundance w/ axionCAMB
for o_idx, o_val in enumerate(omega_ax):
    print("Running RelicFast + axionCAMB for m_ax = ", M_ax_fixed)
    
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
            "m_ax = 1.0e-22", "m_ax = "+f'{M_ax_fixed:.3e}'
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
    
    axion_rho_of_a.append(
        np.loadtxt(outpath+'axion_grhoax_internal.dat')
    )

        
    os.system('mv ./run.ini ./run_fixed_mass_'+str(o_idx)+'.ini')

textvars = []
colors = sns.color_palette("magma", len(omega_ax))

fig, (ax1, ax2) = plt.subplots(2, 1, 
    sharex=True, 
    figsize=(15, 20), 
    gridspec_kw={'height_ratios': [3, 1]}
)
fig.subplots_adjust(hspace=0)
for o_idx, o_val in enumerate(omega_ax): 
    avals = axion_rho_of_a[o_idx][:, 0] 
    rhovals = axion_rho_of_a[o_idx][:, 1] 

    rhointerp = scipy.interpolate.interp1d(avals, rhovals)
    
    a_osc = np.loadtxt(outpath+"axion_aosc.dat")
    rho_osc = rhointerp(a_osc)
    a_dm = np.geomspace(avals[-50], 10**0, 50)
    rho_dm = rho_osc*np.power(a_osc/a_dm, 3)
    #rho_dm = rhovals[-1]*np.power(avals[-1]/a_dm, 3)

    ax1.plot(
        np.log10(avals), 
        rhovals*np.power(10., -10.), 
        label=r'$\Omega_\phi/\Omega_{\rm d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}', 
        linewidth=5.0, 
        color=colors[o_idx]
    )
    ax2.plot(
        np.log10(avals), 
        np.log10(rhovals/omega_cdm_LCDM), 
        label=r'$\Omega_\phi/\Omega_{\rm d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}', 
        linewidth=5.0, 
        color=colors[o_idx]
    )

    if (o_idx==(len(omega_ax)-1)):
        ax1.plot(
            np.linspace(np.log10(avals[0]), np.log10(avals[-1])+0.5, 10), 
            np.array(10*[rhovals[0]])*np.power(10., -10.), 
            color=colors[o_idx], 
            linestyle='dashed',
            linewidth=5.0, 
            label="Dark Energy $\propto a^0$"
        )   
        ax1.plot(
            np.log10(a_dm), 
            rho_dm*np.power(10., -10.), 
            color=colors[o_idx],
            linestyle='dotted',
            linewidth=5.0, 
            label="Dark Matter $\propto a^{-3}$"
        )
        ax2.plot(
            np.linspace(np.log10(avals[0]), np.log10(avals[-1])+0.5, 10), 
            np.log10(10*[rhovals[0]/omega_cdm_LCDM]), 
            color=colors[o_idx], 
            linestyle='dashed',
            linewidth=5.0, 
            label="Dark Energy $\propto a^0$"
        )   
        ax2.plot(
            np.log10(a_dm), 
            np.log10(rho_dm/omega_cdm_LCDM), 
            color=colors[o_idx],
            linestyle='dotted',
            linewidth=5.0, 
            label="Dark Matter $\propto a^{-3}$"
        )
    else:
        pass 
        #ax1.plot(
        #    np.linspace(np.log10(avals[0]), np.log10(avals[-1])+0.5, 10),
        #    np.array(10*[rhovals[0]])*np.power(10., -10.),
        #    color=colors[o_idx],
        #    linestyle='dashed',
        #    linewidth=5.0
        #) 
        ax1.plot(
            np.log10(a_dm), 
            rho_dm*np.power(10., -10.), 
            color=colors[o_idx],
            linewidth=5.0, 
            linestyle='dotted'
        )
        #ax2.plot(
        #    np.linspace(np.log10(avals[0]), np.log10(avals[-1])+0.5, 10),
        #    np.log10(10*[rhovals[0]/omega_cdm_LCDM]),
        #    color=colors[o_idx],
        #    linewidth=5.0, 
        #    linestyle='dashed'
        #) 
        ax2.plot(
            np.log10(a_dm), 
            np.log10(rho_dm/omega_cdm_LCDM), 
            color=colors[o_idx],
            linewidth=5.0, 
            linestyle='dotted'
        )


    #textvar = plt.text(10**-3, rhovals[0], r'$\omega_{\chi, 0} = $'+f'{rho_dm[-1]:.3e}')
    #textvars.append(textvar)
#print(a_osc, avals[-1])

ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax1.set_xlabel(r'$\log(a)$', fontsize=30)
ax1.set_ylabel(r'$10^{-10} \times \omega_\phi$', fontsize=30)
ax1.set_xlim((-6, 0.1))
ax1.tick_params(axis='both', labelsize=25)
#ax1.set_title(r'$M_\phi = 10^{-26}$ eV', fontsize=30)
ax1.legend(fontsize=20, loc='upper right')
ax1.grid(False, which='both', axis='both')

ax2.set_xscale('linear')
ax2.set_yscale('linear')
ax2.set_xlabel(r'$\log(a)$', fontsize=30)
ax2.set_ylabel(r'$\log(\omega_\phi/\omega_{\rm d})$', fontsize=30)
ax2.set_xlim((-6, 0.1))
ax2.tick_params(axis='both', labelsize=25)
#ax2.set_title(r'$M_\chi = $'+f'{M_ax_fixed:.3e}', fontsize=15)
#ax2.set_legend(fontsize=15, loc='upper left')
ax2.grid(False, which='both', axis='both')

plt.savefig(rfpath+"plots/Figure_2.png") 


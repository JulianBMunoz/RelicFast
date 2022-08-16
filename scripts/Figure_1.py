import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()
sns.set_style(style='white')

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
acpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/axionCAMB_Current/"
acpath_output = "Boltzmann_2/transfer_files_0/"
outpath = "/Users/nicholasdeporzio/Downloads/"

M_nu = 0. # Units: eV
redshift = 0.7
kmin = 5.0e-5
kmax = 1.5
Nk = 100
M_halo = 1.0e13

h_lcdm = 0.67
omega_nu = M_nu/93.2
omega_cdm_LCDM = 0.12
omega_b_LCDM = 0.022
Omega_M_LCDM = (omega_b_LCDM + omega_cdm_LCDM)/np.power(h_lcdm, 2.)

Omega_ax = (
    np.array([0.01, 0.02])#([1e-9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10])
    *omega_cdm_LCDM
    /np.power(h_lcdm, 2.)
)
omega_ax = Omega_ax*np.power(h_lcdm, 2.)
omega_cdm = omega_cdm_LCDM - omega_ax

f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

M_ax = np.logspace(-32.5, -22.0, 8)#, 105)
axion_aosc = np.zeros((len(M_ax), len(omega_ax)))
axion_kfs = np.zeros((len(M_ax), len(omega_ax)))

# Calculate Pmm using axionCAMB
os.chdir(acpath)
for m_idx, m_val in enumerate(M_ax):
    print("Running axionCAMB for m_ax = ", m_val)
    for o_idx, o_val in enumerate(omega_ax):
        print("\t For omega_ax = ", o_val)

        axion_kfs[m_idx, o_idx] = (
            np.pi 
            * np.sqrt(m_val*1.56*np.power(10., 29)) 
            * np.power((h_lcdm/2997.)*np.power(1.+redshift, 3.), 0.5)
        ) 
        # axionCAMB for each mass/abundance choice
        reading_file = open(acpath+"params_base.ini", "r")
        new_file_content = ""
        for line in reading_file:
            stripped_line = line.strip()
            new_line = str(stripped_line)
            #new_line = stripped_line.replace("mnu1 = 0.0", "mnu1 = "+f'{M_nu:.2e}')
            #if (M_nu!=0.0):
            #    new_line = new_line.replace(
            #        "tag_sterile_nu = 0", "tag_sterile_nu = 1"
            #    )
            new_file_content += "    " + new_line +"\n"
        new_file_content += "    " + "ombh2 = " + f"{omega_b_LCDM:.4f}" +"\n"
        new_file_content += "    " + "omch2 = " + f"{omega_cdm[o_idx]:.4f}" +"\n"
        new_file_content += "    " + "omnuh2 = " + f"{omega_nu:.4f}" +"\n"
        new_file_content += "    " + "hubble = " + f"{100.*h_lcdm:.4f}" +"\n"
        new_file_content += "    " + "omaxh2 = " + f"{o_val:.4e}" +"\n"
        new_file_content += "    " + "m_ax = " + f"{m_val:.4e}" +"\n"
        new_file_content += "    " + "massless_neutrinos = 3.046"+"\n"
        new_file_content += "    " + "nu_mass_eigenstates = 1"+"\n"
        new_file_content += "    " + "massive_neutrinos = 0"+"\n"
        new_file_content += "    " + "scalar_amp(1) = 2.20e-09"+"\n"
        new_file_content += "    " + "scalar_spectral_index(1) = 0.9655"+"\n"
        new_file_content += "    " + "transfer_num_redshifts = 1"+"\n"
        new_file_content += "    " + "transfer_filename(1) = transfer_out_z" + f"{redshift:.3f}" +"\n"
        new_file_content += "    " + "transfer_redshift(1) = " + f"{redshift:.3f}" +"\n"
        new_file_content += "    " + "output_root = Boltzmann_2/transfer_files_0/"+"\n"

        reading_file.close()
        writing_file = open(acpath+"params_collapse_"+f"{m_idx+o_idx:d}"+".ini", "w")
        writing_file.write(new_file_content)
        writing_file.close()

        os.system('./camb '+"params_collapse_"+f"{m_idx+o_idx:d}"+".ini")

        axion_aosc[m_idx, o_idx] = np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_aosc.dat")

axion_zosc = ((1./axion_aosc)-1.)

np.savetxt(rfpath+"plots/Figure_1_m_ax.txt", M_ax)
np.savetxt(rfpath+"plots/Figure_1_omega_ax.txt", omega_ax)
np.savetxt(rfpath+"plots/Figure_1_axion_aosc.txt", axion_aosc)
np.savetxt(rfpath+"plots/Figure_1_axion_zosc.txt", axion_zosc)
np.savetxt(rfpath+"plots/Figure_1_axion_kfs.txt", axion_kfs)

plot_x = np.geomspace(np.min(M_ax), np.max(M_ax), 100) #Units: [h Mpc^-1]
colors = sns.color_palette("magma", len(omega_ax))
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
for o_idx, o_val in enumerate(omega_ax):
    osc_interp = scipy.interpolate.interp1d(np.log10(M_ax), np.log10(axion_zosc[:,o_idx])) 
    plot_y = osc_interp(np.log10(plot_x))
 
    ax.plot(
        plot_x,
        plot_y,
        label=r'$\Omega_{\phi}/\Omega_{d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}',
        color=colors[o_idx],
        linewidth=5.0
    )

ax.set_xscale('log')
ax.set_yscale('log') 
ax.set_xlabel(r'$M_\phi$ [eV]', fontsize=30)
ax.set_ylabel(r'$z_{\rm osc}$', fontsize=30)
ax.grid(False, which='both', axis='both')
ax.tick_params(axis='both', labelsize=25)
ax.legend(fontsize=25)
plt.savefig(rfpath+"/plots/Figure_1.png")


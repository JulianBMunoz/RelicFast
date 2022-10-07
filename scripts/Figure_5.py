######################################################
####         USER INPUTS
######################################################

import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc

sns.set()
sns.set_style(style='white')
rc('font', **{'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
matplotlib.rcParams.update({
    "font.weight" : "bold",
    "font.size" : 110,
    "axes.labelsize" : 110,
    "axes.labelpad" : 8.0,  
    "xtick.labelsize" : 60, 
    "ytick.labelsize" : 60, 
    "legend.fontsize" : 60, 
    "figure.dpi" : 100, 
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True, 
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

    
axioncamb_data_pk = []
axioncamb_data_tf = []

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
acpath = rfpath+"axionCAMB_Current/"
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

Omega_ax = np.array([1e-9, 0.01, 0.04, 0.07])*omega_cdm_LCDM/np.power(h_lcdm, 2.)
omega_ax = Omega_ax*np.power(h_lcdm, 2.)
omega_cdm = omega_cdm_LCDM - omega_ax

f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

M_ax = [1.0e-30, 1.0e-26]

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

tflookup = {
    "k/h" : 0,
    "cdm" : 1,
    "baryon" : 2,
    "photon" : 3,
    "radiation" : 4,
    "massive nu" : 5,
    "axion" : 6
}

# Calculate Pmm using axionCAMB
os.chdir(acpath)
for m_idx, m_val in enumerate(M_ax):
    print("Running axionCAMB for m_ax = ", m_val)
    for o_idx, o_val in enumerate(omega_ax):
        print("\t For omega_ax = ", o_val)

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

        axioncamb_data_pk.append(
            np.loadtxt(acpath+acpath_output+"_matterpower_1.dat", skiprows=1)
        )
        axioncamb_data_tf.append(
            np.loadtxt(acpath+acpath_output+"_transfer_out_z"+f"{redshift:.3f}", skiprows=1)
        )


plot_x = np.geomspace(10**-3, 10**0.3, 100) #Units: [h Mpc^-1]
colors = sns.color_palette("magma", len(omega_ax))
fig, ax = plt.subplots(1, 1)
linestyles=["solid", "dashed"]
for m_idx, m_val in enumerate(M_ax):
    #ref = np.loadtxt(rfpath+"scripts/2104.07802_FIG2_REF_"+str(m_idx)+".csv", delimiter=",")
    for o_idx, o_val in enumerate(omega_ax):
        idx = len(omega_ax)*m_idx + o_idx

        if o_idx==0:
            TF_lcdm = axioncamb_data_tf[idx]
            TF_lcdm_k = TF_lcdm[:, tflookup['k/h']]*h_lcdm # Careful of units 
            TF_lcdm_cdm = TF_lcdm[:, tflookup['cdm']]*np.power(TF_lcdm_k, 2.)
            TF_lcdm_b = TF_lcdm[:, tflookup['baryon']]*np.power(TF_lcdm_k, 2.)
            TF_lcdm_x = TF_lcdm[:, tflookup['axion']]*np.power(TF_lcdm_k, 2.)
            TF_lcdm_cdm_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_cdm)
            TF_lcdm_b_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_b)
            TF_lcdm_x_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_x)
        else:
            TF = axioncamb_data_tf[idx]
            TF_k = TF[:, tflookup['k/h']]*h_lcdm # Careful of units
            TF_cdm = TF[:, tflookup['cdm']]*np.power(TF_k, 2.)
            TF_b = TF[:, tflookup['baryon']]*np.power(TF_k, 2.)
            TF_x = TF[:, tflookup['axion']]*np.power(TF_k, 2.)
            TF_cdm_interp = scipy.interpolate.interp1d(TF_k, TF_cdm)
            TF_b_interp = scipy.interpolate.interp1d(TF_k, TF_b)
            TF_x_interp = scipy.interpolate.interp1d(TF_k, TF_x)

            Pmm_vals = (
                Pprim(plot_x)
                * (
                    np.power(f_b[m_idx], 2.)*np.power(TF_b_interp(plot_x), 2.)
                    + 2.*f_b[m_idx]*f_cdm[m_idx]*TF_b_interp(plot_x)*TF_cdm_interp(plot_x)
                    + np.power(f_cdm[m_idx], 2.)*np.power(TF_cdm_interp(plot_x), 2.)
                )
            )
    
            Pmm_lcdm_vals = (
                Pprim(plot_x)
                * (
                    np.power(f_b[m_idx], 2.)*np.power(TF_lcdm_b_interp(plot_x), 2.)
                    + 2.*f_b[m_idx]*f_cdm[m_idx]*TF_lcdm_b_interp(plot_x)*TF_lcdm_cdm_interp(plot_x)
                    + np.power(f_cdm[m_idx], 2.)*np.power(TF_lcdm_cdm_interp(plot_x), 2.)
                )
            )

            plot_y = Pmm_vals/Pmm_lcdm_vals

            ax.plot(
                plot_x,
                plot_y,
                label=r'$\omega_{\phi, 0}/\omega_{{\rm d}, 0} = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}'+r"$",
                color=colors[o_idx],
                linestyle=linestyles[m_idx], 
                linewidth=5.0
            )

    ax.set_xscale('log')
    ax.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$')
    ax.set_ylabel(r'$P_{\rm m}/P_{\mathrm{m, }\Lambda \mathrm{CDM}}$')
    ax.grid(False, which='both', axis='both')
    ax.tick_params(axis='both')

    if (m_idx==0):
        ax.legend()
#        ax[m_idx].title((r'$M_\chi = 10^{'+f'{np.log10(m_val):.1f}'+r'}$ eV'))

plt.savefig(rfpath+"plots/Figure_5.png")

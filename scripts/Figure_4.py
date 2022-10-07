######################################################
####         USER INPUTS
######################################################

import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
#import subprocess
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
    'figure.subplot.left': 0.13,
    'figure.subplot.right': 0.98,
    'figure.subplot.top': 0.98,
    'figure.subplot.bottom': 0.06,
#    "figure.constrained_layout.use" : True,
#    "figure.subplot.hspace": 0.001,
    "savefig.pad_inches" : 0.1

})

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"
ac_outpath = "/Users/nicholasdeporzio/Downloads/"

m_ax = np.array([np.power(10., -22.), np.power(10., -31.)])
h = 0.70148
omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
omega_ax = 0.05*omega_cdm_LCDM
Nk=50

omega_cdm = omega_cdm_LCDM-omega_ax

omega_m_LCDM = omega_cdm_LCDM + omega_b_LCDM
H0 = 100.*h
Omega_M_LCDM = omega_m_LCDM/np.power(h, 2.)  
Mnu=0.

################################
rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"
a_val = []
aeq = []
aosc = []

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

for ax_idx, ax_val in enumerate(m_ax):
    print("Running RelicFast + axionCAMB for m_ax = ", ax_val, ", z = 0.65")

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
            "mnu1 = 0.0", "mnu1 = 0."
        ).replace(
            "mnu2 = 0.0", "mnu2 = 0."
        ).replace(
            "m_SN = 0.02", "m_SN = 0."
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = 0.65"
        ).replace(
            "z_collapse_top = 1.4", "z_collapse_top = 1.65"
        ).replace(
            "N_zcoll = 1", "N_zcoll = 1"
        ).replace(
            "hubble = 0.701", "hubble = "+f'{h:.6f}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax:.6e}'
        ).replace(
            "N_klong = 1", "N_klong = "+str(Nk)
        ).replace(
            "omegac = 0.11271", "omegac = "+f'{omega_cdm:.6e}'
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

    print("Generating new data...")

    # Re-compile if axionCAMB output lines don't match RelicFast expectation
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

    if lines_match_flag==True:
        aeq.append(float(np.loadtxt(ac_outpath+"axion_aeq.dat")))
        aosc.append(float(np.loadtxt(ac_outpath+"axion_aosc.dat")))    

#####################################
aeq = np.array(aeq)
aosc = np.array(aosc) 
print(aeq)
print(aosc) 

a_vals = np.logspace(-8., 0., 801) 
zeq = (1./aeq)-1.

Omega_r_LCDM = Omega_M_LCDM/zeq   
omega_r_LCDM = Omega_r_LCDM*np.power(h, 2.)  

def H(H0, Og, Om, a):
    Ol = 1.-Om-Og 
    return (H0*np.sqrt(Og*np.power(a, -4.)+Om*np.power(a, -3.)+Ol))

colors = sns.color_palette("magma", 5)
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
fig.subplots_adjust(hspace=0)
xplot = np.array(a_vals) 

ref1 = H(H0, Omega_r_LCDM[0], Omega_M_LCDM, aosc[0])
ref2 = aosc[0]*np.sqrt(H(H0, Omega_r_LCDM[0], Omega_M_LCDM, aosc[0]))
ref3 = aosc[0]*(H(H0, Omega_r_LCDM[0], Omega_M_LCDM, aosc[0]))
ref4 = aosc[0] 
ref5 = np.sqrt(H(H0, Omega_r_LCDM[0], Omega_M_LCDM, aosc[0]))
ref6 = np.power(aosc[0], -1.)
ref7 = 0.5*np.power(aosc[0], -1.)
ref8 = 2.0*np.power(aosc[0], -1.)
yplot1 = H(H0, Omega_r_LCDM[0], Omega_M_LCDM, xplot)/ref1 
yplot2 = xplot*np.sqrt(H(H0, Omega_r_LCDM[0], Omega_M_LCDM, xplot))/ref2 
yplot3 = xplot*H(H0, Omega_r_LCDM[0], Omega_M_LCDM, xplot)/ref3
yplot4 = xplot/ref4
yplot5 = np.sqrt(H(H0, Omega_r_LCDM[0], Omega_M_LCDM, xplot))/ref5
yplot6 = (np.power(xplot, -1.)/ref6)
yplot7 = (np.power(xplot, -1.)/ref7)
yplot8 = (np.power(xplot, -1.)/ref8)
#ax1.plot(xplot, yplot4, label=r"a", color=colors[0])
#ax1.plot(xplot, yplot2, label=r"$a \sqrt{H}$", color=colors[1])
#ax1.plot(xplot, yplot3, label=r"$aH$", color=colors[2])
ax1.plot(xplot, yplot1, linewidth=3., label=r"$H$", color=colors[4]) 
ax1.plot(xplot, yplot5, linewidth=3., label=r"$k_{\rm J} \propto \sqrt{H}$", color=colors[3], zorder=10)
ax1.plot(xplot, yplot6, linewidth=5., label=r"$k_{\rm m} \propto a^{-1}$", color='black', zorder=1) 
ax1.plot(xplot, yplot7, linewidth=3., color='grey', alpha=0.7) 
ax1.plot(xplot, yplot8, linewidth=3., color='grey', alpha=0.7) 
ax1.plot([aeq[0], aeq[0]], [min(yplot1), max(yplot1)], 
    linewidth=5., color='black', linestyle="dashed", label=r"$a_{\rm eq}$")
ax1.plot([aosc[0], aosc[0]], [min(yplot1), max(yplot1)], 
    linewidth=5., color='black', linestyle="dotted", label=r"$a_{\rm osc}$")
ax1.plot([min(xplot), max(xplot)], [1., 1.], color='black', linestyle='dashdot', label=r"$m_\phi$")
ax1.tick_params(axis='both')
ax1.set_xscale('log') 
ax1.set_yscale('log') 
ax1.set_xlim((min(xplot), max(xplot)))
ax1.set_ylim((1.0e-8, 1.0e1))
ax1.set_yticks(np.logspace(-7, 1, 5))


ref1 = H(H0, Omega_r_LCDM[1], Omega_M_LCDM, aosc[1])
ref2 = aosc[1]*np.sqrt(H(H0, Omega_r_LCDM[1], Omega_M_LCDM, aosc[1]))
ref3 = aosc[1]*(H(H0, Omega_r_LCDM[1], Omega_M_LCDM, aosc[1]))
ref4 = aosc[1] 
ref5 = np.sqrt(H(H0, Omega_r_LCDM[1], Omega_M_LCDM, aosc[1]))
ref6 = np.power(aosc[1], -1.)
ref7 = 0.5*np.power(aosc[1], -1.)
ref8 = 2.0*np.power(aosc[1], -1.)
ref9 = 10.0*np.power(aosc[1], -1.)
yplot1 = H(H0, Omega_r_LCDM[1], Omega_M_LCDM, xplot)/ref1 
yplot2 = xplot*np.sqrt(H(H0, Omega_r_LCDM[1], Omega_M_LCDM, xplot))/ref2 
yplot3 = xplot*H(H0, Omega_r_LCDM[1], Omega_M_LCDM, xplot)/ref3
yplot4 = xplot/ref4
yplot5 = np.sqrt(H(H0, Omega_r_LCDM[1], Omega_M_LCDM, xplot))/ref5
yplot6 = (np.power(xplot, -1.)/ref6)
yplot7 = (np.power(xplot, -1.)/ref7)
yplot8 = (np.power(xplot, -1.)/ref8)
yplot9 = (np.power(xplot, -1.)/ref9)
ax2.plot(xplot, yplot1, linewidth=3., label=r"$H$", color=colors[4]) 
ax2.plot(xplot, yplot5, linewidth=3., label=r"$k_{\rm J} \propto \sqrt{H}$", color=colors[3])
ax2.plot(xplot, yplot6, linewidth=5., label=r"$k_{\rm m} \propto a^{-1}$", color='black') 
ax2.plot(xplot, yplot7, linewidth=3., color='grey', alpha=0.7) 
ax2.plot(xplot, yplot8, linewidth=3., color='grey', alpha=0.7) 
ax2.plot(xplot, yplot9, linewidth=3., color='grey', alpha=0.7) 
ax2.plot([aeq[1], aeq[1]], [1.0e-5, 1.0e15], color='black', 
    linestyle="dashed", linewidth=5., label=r"$a_{\rm eq}$")
ax2.plot([aosc[1], aosc[1]], [1.0e-5, 1.0e15], color='black', 
    linestyle="dotted", linewidth=5., label=r"$a_{\rm osc}$")
ax2.plot([min(xplot), max(xplot)], [1., 1.], color='black', linestyle='dashdot', label=r"$m_\phi$")
#ax2.plot(xplot[0:-1], np.diff(np.log10(yplot1))/np.diff(np.log10(xplot))) 
ax2.tick_params(axis='both')
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlim((min(xplot), max(xplot)))
ax2.set_ylim((1.0e-3, 1.0e6))
ax2.set_yticks(np.logspace(-1, 5, 4))
ax2.legend(loc="lower left")

#fig.supxlabel(r"$a$") 
#fig.supylabel(r"$k_{{\rm physical}} {\rm ~[} m_\phi {\rm ]}$")
#fig.tight_layout()
fig.text(0.001, 0.5, r"$k_{{\rm physical}} {\rm ~[} m_\phi {\rm ]}$", 
    rotation='vertical', va='center')
fig.text(0.555, 0.001, r"$a$", 
    ha='center', rotation='horizontal')

plt.savefig(rfpath+"plots/Figure_4.png")

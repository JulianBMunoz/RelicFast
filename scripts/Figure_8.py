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
use_existing_data=True
data_save_level=2
######################################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"
ac_outpath = "/Users/nicholasdeporzio/Downloads/"

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226

m_ax = np.array([
    #np.power(10., -23.0),
    np.power(10., -25.0),
    np.power(10., -27.0),
    np.power(10., -29.0),
    #np.power(10., -30.0),
    np.power(10., -31.0),
    np.power(10., -32.0),
    #np.power(10., -32.5)
])
omega_ax = 0.05*omega_cdm_LCDM
sum_massive_nu = 0.
redshifts = np.concatenate((np.array([0.01]), np.linspace(1., 10., 10))) 
kmin = 0.9e-4
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

kfs = np.array([[
    np.pi * np.sqrt(m*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+z, 3.), 0.5)
    for z in redshifts]
    for m in m_ax]
) 

######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

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
for ax_idx, ax_val in enumerate(m_ax): 
    for z_idx, z_val in enumerate(redshifts): 
        print("Running RelicFast + axionCAMB for m_ax = ", ax_val, ", z = ", z_val)
        
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
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{z_val:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                #"hubble = 0.701", "hubble = "+f'{h[ax_idx]:.6f}'
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

        # Check if data already exists
        b1e_path = (rfpath+"plots/Figure_8_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt") 
        b1l_path = (rfpath+"plots/Figure_8_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt") 
        ax_background_path = (rfpath+"plots/Figure_8_axbackground_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt")
        ax_aosc_path = (rfpath+"plots/Figure_8_axaosc_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt") 
 
        if (
            os.path.exists(b1e_path)    
            and os.path.exists(b1l_path)
            and os.path.exists(ax_background_path)
            and os.path.exists(ax_aosc_path) 
            and (use_existing_data==True)
        ):
            print("Loading existing data...")
            test_in_b1e = np.loadtxt(b1e_path)
            test_in_b1l = np.loadtxt(b1l_path)
            test_in_ax_background = np.loadtxt(ax_background_path)
            test_in_ax_aosc = float(np.loadtxt(ax_aosc_path))
            print("Loading ", b1e_path)
            print("Loading ", b1l_path)
            print("Loading ", ax_background_path)
            print("Loading ", ax_aosc_path)
            data_b1e=np.array(test_in_b1e)
            data_b1l=np.array(test_in_b1l)
            ax_background = np.array(test_in_ax_background)
            ax_aosc = test_in_ax_aosc
            
        else:
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
    
                    expected_axioncamb_output_lines = int(ac_lines_out)

            if lines_match_flag==True:
                data_b1e = np.loadtxt(
                    'output/result-0/bias_Euler_z'+f'{z_val:.2f}'
                    +'_M13.00_Nk'+str(Nk)
                    +'.dat', skiprows=1
                )
                data_b1l = np.loadtxt(
                    'output/result-0/bias_Lagrangian_z'+f'{z_val:.2f}'
                    +'_M13.00_Nk'+str(Nk)
                    +'.dat', skiprows=1
                )
                ax_background = np.loadtxt(ac_outpath+"axion_background.dat")
                ax_aosc = float(np.loadtxt(ac_outpath+"axion_aosc.dat"))
                print("a_osc = ", ax_aosc)
                if data_save_level>1:
                    np.savetxt((rfpath+"plots/Figure_8_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt"), data_b1e)
                    np.savetxt((rfpath+"plots/Figure_8_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt"), data_b1l) 
                    np.savetxt((rfpath+"plots/Figure_8_axbackground_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt"), ax_background)
                    np.savetxt((rfpath+"plots/Figure_8_axaosc_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{omega_ax:.6f}"+"_z"+f"{z_val:.2f}"+".txt"), [ax_aosc])        

        axion_background.append(ax_background)
        axion_aosc.append(ax_aosc)
        data_eulbias.append(data_b1e)
        data_lagbias.append(data_b1l) 
    
        os.system('mv ./run.ini ./run_'+str(ax_idx)+'.ini')
        os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini')  
kmin = max([min(data_lagbias[idx][:,0]) for idx, val in enumerate(data_lagbias)])
kmax = min([max(data_lagbias[idx][:,0]) for idx, val in enumerate(data_lagbias)])

k_interp = np.geomspace(kmin, kmax, Nk)
k_interp_table = ([[
    scipy.interpolate.interp1d(data_lagbias[m_idx*len(redshifts)+z_idx][:,0], data_lagbias[m_idx*len(redshifts)+z_idx][:,1])(k_interp)
    for z_idx, z_val in enumerate(redshifts)]
    for m_idx, m_val in enumerate(m_ax)]
)
np.shape(k_interp_table)

k_z_interp_table = [
    scipy.interpolate.interp2d(np.log10(k_interp), redshifts, k_interp_table[m_idx])
    for m_idx, m_val in enumerate(m_ax)
]

print(axion_aosc) 
axion_aosc = np.reshape(axion_aosc, (len(m_ax), len(redshifts)))
axion_zosc = (1./axion_aosc)-1.
print(axion_zosc)
colors = sns.color_palette("icefire", 2*len(redshifts)+1)

fig, ax = plt.subplots(len(m_ax), 1,
    sharex=True,
    figsize=(15., 7.5*len(m_ax)),
    gridspec_kw={'height_ratios': [1]*len(m_ax)}
)
fig.subplots_adjust(hspace=0)
for m_idx, m_val in enumerate(m_ax): 
    zosc = axion_zosc[m_idx, z_idx]
    Nzosc = len(redshifts[redshifts<zosc])
    Nznonosc = len(redshifts[redshifts<zosc])
    mididx = len(redshifts)
    startidx = mididx + Nzosc
    
    for z_idx, z_val in enumerate(redshifts): 
        norm = (k_z_interp_table[m_idx](np.log10(kmin), z_val))
        ax[m_idx].plot(
            k_interp, 
            k_z_interp_table[m_idx](np.log10(k_interp), z_val)/norm, 
            color=colors[startidx-z_idx], 
            label=r"$z = "+f"{z_val:.2f}"+r"$", 
            zorder=(m_idx*len(redshifts)+z_idx),
            alpha=0.6, 
            linewidth=2.
        )
        pt = ax[m_idx].scatter(
            kfs[m_idx][z_idx], 
            k_z_interp_table[m_idx](np.log10(kfs[m_idx][z_idx]), z_val)/norm, 
            facecolors='none',
            edgecolors=colors[startidx-z_idx],
            marker="o",
            linewidth=2., 
            s=200, 
            zorder=(len(m_ax)*len(redshifts)+5)
        )
        ax[m_idx].set_xscale('log')
        ax[m_idx].tick_params(axis='y', labelsize=30)
        ax[m_idx].text(
            np.power(10., -1.3), 1., 
            r"$m_\phi = 10^{"+f"{np.log10(m_val):.1f}"+r"}$ eV", 
            fontsize=30, 
            bbox=dict(facecolor='white', edgecolor='black', pad=10.0)
        ) 
        if (z_idx==(len(redshifts)-1)):
            ax[m_idx].plot(
                [0.70148*0.015, 0.70148*0.015],
                [0.9, 1.1],
                color='red',
                label=r'$k_{\rm eq}$',
                linewidth=2.
            )
        ax[m_idx].set_xlim((kmin, kmax))
        ax[m_idx].set_ylim((0.997, 1.04))
        if ((m_idx==0) and (z_idx==0)):
            pass 
            #pt.set_label(r"$k_{\rm fs}$")
        #if m_idx==0:
        #    ax[m_idx].legend(fontsize=30, loc="upper left")
        if (m_idx==(len(m_ax)-1)): 
            ax[m_idx].set_xlabel(r"k [Mpc$^{-1}$]", fontsize=40)
        if (m_idx==int(len(m_ax)/2)): 
            ax[m_idx].set_ylabel(r"$b_1^L(k)~/~b_1^L(k_{*})$", fontsize=40)            
fig.subplots_adjust(hspace=0)
plt.xticks(fontsize=30)
plt.savefig(rfpath+"plots/Figure_8.png")

colors = sns.color_palette("icefire", 2*len(redshifts)+1)

fig, ax = plt.subplots(len(m_ax), 1,
    sharex=True,
    figsize=(15., 7.5*len(m_ax)),
    gridspec_kw={'height_ratios': [1]*len(m_ax)}
)
fig.subplots_adjust(hspace=0)
for m_idx, m_val in enumerate(m_ax): 
    zosc = axion_zosc[m_idx, z_idx]
    Nzosc = len(redshifts[redshifts<zosc])
    Nznonosc = len(redshifts[redshifts<zosc])
    mididx = len(redshifts)
    startidx = mididx + Nzosc
    
    for z_idx, z_val in enumerate(redshifts): 
        #norm = (k_z_interp_table[m_idx](np.log10(kmin), z_val))
        norm = 1. 
        ax[m_idx].plot(
            k_interp, 
            k_z_interp_table[m_idx](np.log10(k_interp), z_val)/norm, 
            color=colors[startidx-z_idx], 
            label=r"$z = "+f"{z_val:.2f}"+r"$", 
            zorder=(m_idx*len(redshifts)+z_idx),
            alpha=0.6
        )
        pt = ax[m_idx].scatter(
            kfs[m_idx][z_idx], 
            k_z_interp_table[m_idx](np.log10(kfs[m_idx][z_idx]), z_val)/norm, 
            marker=".", 
            s=50, 
            color=colors[startidx-z_idx],
            zorder=(len(m_ax)*len(redshifts)+5)
        )
        ax[m_idx].set_xscale('log')
        ax[m_idx].set_yscale('log')
        ax[m_idx].grid(True, which='major')
        #ax[m_idx].text(
        #    np.power(10., -0.6), 1., 
        #    r"$m_\phi = 10^{"+f"{np.log10(m_val):.0f}"+r"}$ eV", 
        #    fontsize=15, 
        #    bbox=dict(facecolor='white', edgecolor='black', pad=10.0)
        #)  
        #if (z_idx==(len(redshifts)-1)):
        #    ax[m_idx].plot(
        #        [0.70148*0.015, 0.70148*0.015],
        #        [1., 1.035],
        #        color='red',
        #        label=r'$k_{\rm eq}$',
        #        linewidth=1.
        #    )
        ax[m_idx].set_xlim((kmin, 1.4*kmax))
        if ((m_idx==0) and (z_idx==0)): 
            pt.set_label(r"$k_{\rm fs}$")
        if m_idx==0:
            ax[m_idx].legend(fontsize=8, loc="upper left")
        if (m_idx==(len(m_ax)-1)): 
            ax[m_idx].set_xlabel(r"k [Mpc$^{-1}$]", fontsize=15)
        if (m_idx==int(len(m_ax)/2)): 
            ax[m_idx].set_ylabel(r"$b_1^L(k)$", fontsize=15)
            
fig.subplots_adjust(hspace=0)
plt.savefig(rfpath+"plots/Figure_8b.png")



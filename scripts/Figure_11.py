# NOTES! Output will be save in the `plots/` directory of your
# RelicFast install location, which you should specify below.
 
######################################################
####    EDIT THESE LINES FOR YOUR SYSTEM 
######################################################

#Path to RelicFast install
rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"


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
    "font.size" : 60,
    "axes.labelsize" : 110,
    "axes.labelpad" : 8.0,  
    "xtick.labelsize" : 60, 
    "ytick.labelsize" : 60, 
    "legend.fontsize" : 60, 
    "figure.dpi" : 300, 
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True, 
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})
use_existing_data=True
data_save_level=2

# Set cosmological parameters 
omega_cdm_LCDM = 0.1127 # Units: none
omega_b_LCDM = 0.02226 # Units: none
h_LCDM = 0.70148 # Units: none 
m_ax = np.logspace(-32., -22., 41) # Units: eV
omega_ax = np.concatenate(( # Units: none
    np.array([1.0e-9*omega_cdm_LCDM]), 
    np.linspace(0.010, 0.100, 37)*omega_cdm_LCDM 
))
print("Axion masses: ", m_ax)
print("Axion abundances (% CDM): ", omega_ax/omega_cdm_LCDM)
sum_massive_nu = 0. # Units: eV
omega_nu = sum_massive_nu/93.2 # Units: none
redshift = 0.65 # Units: none
kmin = 1.0e-4 # Units: 
kmax = 1.0 # Units: 
Nk = 50 # Units: none

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

# Scales to compute bias values at 
krefs = np.logspace(-3.5, -0.5, 4) # Units: 

# Scale to normalize step plots to 
knorm = np.power(10., -3.9) # Units:  

# Initialize results arrays
b1e = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1l = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1e_step = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1l_step = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
spontaneous_failures=[]

# Compile RelicFast for use with axions
rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"
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

# Run RelicFast for each axion mass/abundance  
for ax_idx, ax_val in enumerate(m_ax):
    print("Running RelicFast + axionCAMB for m_ax = "+f"{ax_val:.2e}"+"...") 
    for o_idx, o_val in enumerate(omega_ax):
        print("\t for omega_ax = "+f"{o_val:.3f}"+"...")

        # Clear old data 
        os.system('rm -r output/result-0/')
        os.system('rm -r '+'Boltzmann_0')
        os.system('rm -r '+'Boltzmann_1')
        os.system('rm -r '+'Boltzmann_2')
        os.system('rm ./run.ini')

        # Calculate derived quantities
        omega_cdm = omega_cdm_LCDM - o_val
        f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
        f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)
        
        h = h_LCDM*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+o_val)/(omega_b_LCDM+omega_cdm_LCDM))
        kfs = np.pi * np.sqrt(ax_val*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)

        # Make a RelicFast .ini file
        reading_file = open("2006.09395.ini", "r")
        new_file_content = ""
        for line in reading_file:
            stripped_line = line.strip()
            new_line = stripped_line.replace(
                "mnu1 = 0.0", "mnu1 = "+f'{sum_massive_nu:.2e}'
            ).replace(
                "mnu2 = 0.0", "mnu2 = "+f'{sum_massive_nu:.2e}'
            ).replace(
                "m_SN = 0.02", "m_SN = "+f'{sum_massive_nu:.2e}'
            ).replace(
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{o_val:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm:.6e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{ax_val:.6e}'
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

        # Check if data already exists
        b1e_path = (rfpath+"plots/Figure_11_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
            +"_omegaaxion"+f"{o_val:.6f}"+".txt")
        b1l_path = (rfpath+"plots/Figure_11_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
            +"_omegaaxion"+f"{o_val:.6f}"+".txt")

        if (os.path.exists(b1e_path) and os.path.exists(b1l_path) and (use_existing_data==True)): 
            print("Loading existing data...")
            test_in_b1e = np.loadtxt(b1e_path)
            test_in_b1l = np.loadtxt(b1l_path) 
            print("Loading ", b1e_path)
            print("Loading ", b1l_path) 
            data_eulbias=np.array(test_in_b1e)
            data_lagbias=np.array(test_in_b1l)  

        # Generate data if it doesn't already exist 
        else:  
            print("Generating new data...") 
            
            lines_match_flag=False
            while lines_match_flag==False:    
                try:      
                    os.system('./relicfast run.ini')
                    with open(rfpath+"/Boltzmann_2/transfer_files_0/_transfer_out_z200.000", 'r') as fp:
                        ac_lines_out = sum(1 for line in fp)
                    fp.close()
                except: 
                    print("Spontaneous axionCAMB fault!!!")
                    spontaneous_failures.append(int(ax_idx*len(omega_ax)+o_idx)) 
                    break 
                
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
                data_eulbias = np.loadtxt(
                    'output/result-0/bias_Euler_z'+f'{redshift:.2f}'
                    +'_M13.00_Nk'+str(Nk)
                    +'.dat', skiprows=1
                )
                data_lagbias = np.loadtxt(
                    'output/result-0/bias_Lagrangian_z'+f'{redshift:.2f}'
                    +'_M13.00_Nk'+str(Nk)
                    +'.dat', skiprows=1
                )
    
                if data_save_level>1: 
                    np.savetxt((rfpath+"plots/Figure_11_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{o_val:.6f}"+".txt"), data_eulbias)
                    np.savetxt((rfpath+"plots/Figure_11_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                        +"_omegaaxion"+f"{o_val:.6f}"+".txt"), data_lagbias)
            else: 
                data_eulbias = np.array([[1.0e-4, 0.], [1.0, 0.]])
                data_lagbias = np.array([[1.0e-4, 0.], [1.0, 0.]])

        eulbias_interp = scipy.interpolate.interp1d(data_eulbias[:,0], data_eulbias[:,1])
        lagbias_interp = scipy.interpolate.interp1d(data_lagbias[:,0], data_lagbias[:,1]) 

        for kidx, kval in enumerate(krefs): 
            b1e[kidx, ax_idx, o_idx] = eulbias_interp(kval) 
            b1l[kidx, ax_idx, o_idx] = lagbias_interp(kval)
            b1e_step[kidx, ax_idx, o_idx] = eulbias_interp(kval)/eulbias_interp(knorm) 
            b1l_step[kidx, ax_idx, o_idx] = lagbias_interp(kval)/lagbias_interp(knorm) 

        os.system('mv ./run.ini ./run_'+str(ax_idx*len(omega_ax)+o_idx)+'.ini')
        os.system(
            'mv ./axionCAMB_Current/params_collapse.ini'
            +' ./axionCAMB_Current/params_collapse_'+str(ax_idx*len(omega_ax)+o_idx)+'.ini'
        )  

np.savetxt(rfpath+"plots/Figure_11_Failures.txt", np.array(spontaneous_failures, dtype='int'))
print("Spontaneous failures: ", np.array(spontaneous_failures, dtype='int'))

def fmt(x):
    s = f"{(x-1.)*100.:.0f}"
    return rf"${s} \%$" if plt.rcParams["text.usetex"] else f"{s} %"

for kidx, kval in enumerate(krefs):
    if data_save_level>0:  
        np.savetxt(rfpath+"plots/Figure_11_b1e_logk"+f"{np.log10(kval):.3f}"+".txt", b1e[kidx])
        np.savetxt(rfpath+"plots/Figure_11_b1l_logk"+f"{np.log10(kval):.3f}"+".txt", b1l[kidx])
        np.savetxt(rfpath+"plots/Figure_11_b1estep_logk"+f"{np.log10(kval):.3f}"+".txt", b1e_step[kidx])
        np.savetxt(rfpath+"plots/Figure_11_b1lstep_logk"+f"{np.log10(kval):.3f}"+".txt", b1l_step[kidx])

    Z = np.transpose(np.nan_to_num(b1l[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
    fig, ax = plt.subplots(1,1)
    hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto")
    ax.tick_params(axis='both')
    ax.set_xlabel(r"$\log{M_\phi / {\rm [eV]}}$")
    ax.set_ylabel("$\omega_\phi / \omega_d$")
    ax.set_ylim((0.,0.1)) 
    cbar = plt.colorbar(hmap)
    plt.savefig(rfpath+"plots/Figure_11_b1l_logk"+f"{np.log10(kval):.3f}"+".png")

    Z = np.transpose(np.nan_to_num(b1l_step[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
    fig, ax = plt.subplots(1,1)
    hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(b1l_step), vmax=np.max(b1l_step), shading="auto")
    if (np.max(Z)>1.01):
        CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
        ax.clabel(CS, CS.levels, inline=True, fmt=fmt)
    ax.tick_params(axis='both')
    ax.set_xlabel(r"$\log{\left(M_\phi ~/~ {\rm [eV]}\right)}$")
    ax.set_ylabel("$\omega_\phi ~/~ \omega_d$")
    ax.set_ylim((0.,0.1)) 
    cbar = plt.colorbar(hmap)
    cbar.set_label(label=(
        r"$b^1_L(k)~/~b^1_L(k_{\rm ref})$")#, size=50
    )
    plt.savefig(rfpath+"plots/Figure_11_b1lstep_logk"+f"{np.log10(kval):.3f}"+".png")
    if (kidx==(len(krefs)-1)): 
        plt.savefig(rfpath+"plots/Figure_11.png") 

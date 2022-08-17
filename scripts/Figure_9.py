# NOTES! Output will be save in the `plots/` directory of your
# RelicFast install location, which you should specify below.
 
######################################################
####    EDIT THESE LINES FOR YOUR SYSTEM 
######################################################

#Path to RelicFast install
rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"


######################################################
####    INTERNAL - DO NOT EDIT BELOW HERE 
######################################################

# Import existing data or generate new data
use_existing_data = True 
data_save_level = 2

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D

# Plot settings
sns.set()
sns.set_style("white")

# Set cosmological parameters 
omega_cdm_LCDM = 0.1127 # Units: none
omega_b_LCDM = 0.02226 # Units: none
h_LCDM = 0.70148 # Units: none 
m_ax = np.logspace(-32., -22., 21) # Units: eV
omega_ax = np.concatenate(( # Units: none
    np.array([1.0e-9*omega_cdm_LCDM]), 
    np.linspace(0.010, 0.100, 19)*omega_cdm_LCDM 
))
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

# Initialize results arrays
b1l = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1e = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none

#if use_existing_data==True: 
#    print("Attempting to load existing data...") 
#    for kidx, kval in enumerate(krefs): 
#        b1e_path = rfpath+"plots/Figure_9_b1e_logk"+f"{np.log10(kval):.3f}"+".txt"
#        b1l_path = rfpath+"plots/Figure_9_b1l_logk"+f"{np.log10(kval):.3f}"+".txt"
#        if (os.path.exists(b1e_path) and os.path.exists(b1l_path)): 
#            test_in_b1e = np.loadtxt(b1e_path)
#            test_in_b1l = np.loadtxt(b1l_path) 
#            if (
#                (np.shape(test_in_b1e)==np.shape(b1e[kidx]))
#                and (np.shape(test_in_b1l)==np.shape(b1l[kidx]))
#            ):
#                print("Loading ", b1e_path)
#                print("Loading ", b1l_path) 
#                b1e[kidx]=np.array(test_in_b1e)
#                b1l[kidx]=np.array(test_in_b1l)   

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

        b1e_path = (rfpath+"plots/Figure_9_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
            +"_omegaaxion"+f"{o_val:.6f}"+".txt")
        b1l_path = (rfpath+"plots/Figure_9_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
            +"_omegaaxion"+f"{o_val:.6f}"+".txt")
        if (os.path.exists(b1e_path) and os.path.exists(b1l_path) and (use_existing_data==True)): 
            print("Loading existing data...")
            test_in_b1e = np.loadtxt(b1e_path)
            test_in_b1l = np.loadtxt(b1l_path) 
            print("Loading ", b1e_path)
            print("Loading ", b1l_path) 
            data_eulbias=np.array(test_in_b1e)
            data_lagbias=np.array(test_in_b1l)  

        else:  
            print("Generating new data...") 
            omega_cdm = omega_cdm_LCDM - o_val
            f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
            f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)
            
            h = h_LCDM*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+o_val)/(omega_b_LCDM+omega_cdm_LCDM))
            kfs = np.pi * np.sqrt(ax_val*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)
            
            # Clear old data 
            os.system('rm -r output/result-0/')
            os.system('rm -r '+'Boltzmann_0')
            os.system('rm -r '+'Boltzmann_1')
            os.system('rm -r '+'Boltzmann_2')
            os.system('rm ./run.ini')
            
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
                np.savetxt((rfpath+"plots/Figure_9_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_omegaaxion"+f"{o_val:.6f}"+".txt"), data_eulbias)
                np.savetxt((rfpath+"plots/Figure_9_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_omegaaxion"+f"{o_val:.6f}"+".txt"), data_lagbias)

        eulbias_interp = scipy.interpolate.interp1d(data_eulbias[:,0], data_eulbias[:,1])
        lagbias_interp = scipy.interpolate.interp1d(data_lagbias[:,0], data_lagbias[:,1]) 

        for kidx, kval in enumerate(krefs): 
            b1e[kidx, ax_idx, o_idx] = eulbias_interp(kval) 
            b1l[kidx, ax_idx, o_idx] = lagbias_interp(kval) 

        os.system('mv ./run.ini ./run_'+str(ax_idx+o_idx)+'.ini')
        os.system(
            'mv ./axionCAMB_Current/params_collapse.ini'
            +' ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini'
        )  

for kidx, kval in enumerate(krefs):
    if data_save_level>0:  
        np.savetxt(rfpath+"plots/Figure_9_b1e_logk"+f"{np.log10(kval):.3f}"+".txt", b1e[kidx])
        np.savetxt(rfpath+"plots/Figure_9_b1l_logk"+f"{np.log10(kval):.3f}"+".txt", b1l[kidx])

    Z = np.nan_to_num(b1l[kidx])
    X, Y = np.meshgrid(m_ax, omega_ax/omega_cdm_LCDM)
    fig, ax = plt.subplots(1,1, figsize=(15, 15))
    hmap = ax.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z), vmax=np.max(Z))
    ax.set_xscale('log')
    ax.tick_params(axis='both', labelsize=25)
    ax.set_xlabel("$M_\phi$ [eV]", fontsize=30)
    ax.set_ylabel("$\Omega_\phi / \Omega_d$", fontsize=30)
    cbar = plt.colorbar(hmap)
    cbar.set_label(label=(r"$b^1_L(k=10^{"+f"{np.log10(kval):.1f}"+r"})$"), size=30)
    cbar.ax.tick_params(labelsize=25) 
    plt.savefig(rfpath+"plots/Figure_9_b1l_logk"+f"{np.log10(kval):.3f}"+".png")

# Can use circle patches if you want to control the width of the plot circles carefully
#https://stackoverflow.com/questions/48172928/scale-matplotlib-pyplot-axes-scatter-markersize-by-x-scale/48174228#48174228

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
    "figure.figsize" : [30, 40],
    'figure.subplot.left': 0.13,
    'figure.subplot.right': 0.98,
    'figure.subplot.top': 0.98,
    'figure.subplot.bottom': 0.06,
    #"figure.constrained_layout.use" : True, 
    #"figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})
use_existing_data=True
data_save_level=2

######################################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

m_ax = np.array([
    10**-24, 
    10**-28, 
    10**-32
])
m_halo = np.logspace(13., 15., 3)

omega_ax = np.array(len(m_ax)*[0.1*omega_cdm_LCDM])

sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-5
kmax = 0.5
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

h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

kfs = np.pi * np.sqrt(m_ax*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)

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
for ax_idx, ax_val in enumerate(m_ax): 
    for mh_idx, mh_val in enumerate(m_halo): 
        print("Running RelicFast + axionCAMB for m_ax = ", ax_val)
        
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
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h[ax_idx]:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax[ax_idx]:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "kbot = 1e-4", "kbot = "+f"{kmin:.2e}"
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm[ax_idx]:.6e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{ax_val:.6e}'
            ).replace(
                "Mhalo_min = 1e13", "Mhalo_min = "+f'{mh_val:.6e}'
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
        b1e_path = (rfpath+"plots/Figure_9_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_mhalo"+f"{mh_val:.2e}"+".txt")
        b1l_path = (rfpath+"plots/Figure_9_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_mhalo"+f"{mh_val:.2e}"+".txt")

        if (os.path.exists(b1e_path) and os.path.exists(b1l_path) and (use_existing_data==True)): 
            print("Loading existing data...")
            test_in_b1e = np.loadtxt(b1e_path)
            test_in_b1l = np.loadtxt(b1l_path) 
            print("Loading ", b1e_path)
            print("Loading ", b1l_path) 
            data_eulbias.append(test_in_b1e)
            data_lagbias.append(test_in_b1l)  

        else:  
            print("Generating new data...") 
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
            #pm_idx = np.argmin(np.abs(z_vals - 0.0))
            print('Requested/found redshift: ', redshift, z_vals[pm_idx])
            
            print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
            
            #data_pm.append(
            #    np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            #)
            #
            #data_tf.append(
            #    np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
            #)
            
            data_eulbias.append(
                np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_vals[pm_idx]:.2f}'+'_M'+f'{np.log10(mh_val):.2f}'+'_Nk'+str(Nk)+'.dat', skiprows=1)
            )
            data_lagbias.append(
                np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_vals[pm_idx]:.2f}'+'_M'+f'{np.log10(mh_val):.2f}'+'_Nk'+str(Nk)+'.dat', skiprows=1)
            )
    
            if data_save_level>1:
                np.savetxt((rfpath+"plots/Figure_9_b1e_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_mhalo"+f"{mh_val:.2e}"+".txt"), data_eulbias[-1])
                np.savetxt((rfpath+"plots/Figure_9_b1l_logmaxion"+f"{np.log10(ax_val):.3f}"
                    +"_mhalo"+f"{mh_val:.2e}"+".txt"), data_lagbias[-1])
    
            os.system('mv ./run.ini ./run_'+str(mh_idx+len(m_ax)*ax_idx)+'.ini')
            os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(mh_idx+len(m_ax)*ax_idx)+'.ini')  



# eulerian bias plots
#plt.figure()
#plt.xscale('log')
#for ax_idx, ax_val in enumerate(m_ax): 
#    if ((ax_idx==0) or (ax_idx==(len(m_ax)-1))): 
#        continue 
#    
#    kvals = data_eulbias[ax_idx][:, 0]
#    eulbiasvals = data_eulbias[ax_idx][:, 1]
#    
#    eulbiasinterp = scipy.interpolate.interp1d(kvals, eulbiasvals)
#    eulbiasplot = eulbiasinterp(kplot)
#    
#    rgb = np.zeros(3)
#    rgb[0] = 3./255.
#    rgb[1] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
#    rgb[2] = (len(m_ax)-1-ax_idx) * (255/(len(m_ax)-1)) / 255.
#    
#    yplot = eulbiasplot/eulbiasplot[0]
#    
#    plt.plot(kplot, yplot, label=r'$m_{\chi, i}=$'+f'{ax_val:.3e}'+r' eV', color=tuple(rgb))
#    plt.plot([kfs[ax_idx], kfs[ax_idx]], [min(yplot), max(yplot)], color=tuple(rgb), linestyle='dashed')
#    
#plt.plot([0.7*0.015, 0.7*0.015], [1., 1.01], color='red', label=r'$k_{eq}$')
#plt.xlabel(r'$k ~[{\rm Mpc}^{-1}]$')
#plt.ylabel(r'$b_1(k)/b_1(k_{\rm ref})$')
#plt.title(r'$\omega_\chi = 0.05\times\omega_{cdm},  ~z = 0.65, ~\Sigma M_\nu = 0$ eV')
##plt.legend(lines, labels)
#plt.legend()
##plt.yscale('log')
#plt.grid(False, which='both', axis='both')

# bias plots 
kplot = np.geomspace(10**-4.0, 0.6, 100)
lstyles = ["solid", "dashed", "dotted"]
colors=sns.color_palette("magma", len(m_ax))
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
fig.subplots_adjust(hspace=0)
for ax_idx, ax_val in enumerate(m_ax): 
    for mh_idx, mh_val in enumerate(m_halo):
        didx = (mh_idx+ax_idx*len(m_ax)) 
        lagkvals = data_lagbias[didx][:, 0]
        eulkvals = data_eulbias[didx][:, 0]

        lagbiasvals = data_lagbias[didx][:, 1]
        eulbiasvals = data_eulbias[didx][:, 1]
        
        lagbiasinterp = scipy.interpolate.interp1d(lagkvals, lagbiasvals)
        eulbiasinterp = scipy.interpolate.interp1d(eulkvals, eulbiasvals)

        lagbiasplot = lagbiasinterp(kplot)
        eulbiasplot = eulbiasinterp(kplot)
        
        yplot1 = lagbiasplot/lagbiasplot[0]
        yplot2 = eulbiasplot/eulbiasplot[0]

#        aeq = 3.0988218673675137E-004
#        mosc = np.loadtxt(rfpath+"plots/Figure_1_m_ax.txt")
#        aosc = np.loadtxt(rfpath+"plots/Figure_1_axion_aosc.txt")[:,0]
#        logaosc_interp = scipy.interpolate.interp1d(np.log10(mosc), np.log10(aosc))
#        aoscval = np.power(10., logaosc_interp(np.log10(ax_val)))
#        RLCDM = 1. + 0.0048*np.tanh(4.*kplot/(0.70148*0.015)) #LCDM bias step 
#        RPHI = 1.+1.*np.tanh(np.power((aeq/aoscval)*(1.-aoscval), 1.))*(omega_ax[ax_idx]/(omega_cdm_LCDM+omega_b_LCDM))*np.tanh(kplot/kfs[ax_idx])  
#        yplot3 = RLCDM*RPHI        

        if mh_idx==0:         
            ax1.plot(
                kplot, 
                yplot1, 
                label=r'$m_{\phi}= 10^{'+f'{np.log10(ax_val):.0f}'+r'} {\rm ~eV}$', 
                color=colors[ax_idx], 
                linewidth=5.,
                linestyle=lstyles[mh_idx])
        if ax_idx==0: 
            ax2.plot(
                kplot, 
                yplot2, 
                color=colors[ax_idx], 
                linewidth=5.,
                label=r"$M= 10^{"+f"{np.log10(mh_val):.0f}"+r"}~M_\odot$", 
                linestyle=lstyles[mh_idx])
        ax1.plot(
            kplot, 
            yplot1, 
            color=colors[ax_idx], 
            linewidth=5.,
            linestyle=lstyles[mh_idx])
        ax2.plot(
            kplot, 
            yplot2, 
            color=colors[ax_idx], 
            linewidth=5.,
            linestyle=lstyles[mh_idx])
#        ax1.plot(
#            kplot, 
#            yplot3, 
#            color=colors[ax_idx], 
#            linewidth=3., 
#            linestyle='dotted'
#        ) 
        #try: 
        #    ax.scatter(
        #        kfs[ax_idx],
        #        lagbiasinterp(kfs[ax_idx])/lagbiasplot[0], 
        #        facecolors='none', 
        #        edgecolors=tuple(rgb),
        #        marker='o', 
        #        linewidth=3., 
        #        s=400, 
        #        zorder=3)
        #except:
        #    pass 
        
    #ax.plot(
    #    [0.70148*0.015, 0.70148*0.015], 
    #    [1., 1.07], 
    #    color='red', 
    #    label=r'$k_{\rm eq}$',
    #    linewidth=5.)

ax1.set_xscale('log')
ax1.set_ylabel(r'$b_1^L(k)/b_1^L(k_{\rm ref})$')
ax1.tick_params(axis='both')
ax1.set_xlim((1.0e-4, max(kplot)))
ax1.set_ylim((0.997, 1.059))
ax1.legend()
ax1.grid(False, which='both', axis='both')
ax2.set_xscale('log')
ax2.set_ylabel(r'$b_1(k)/b_1(k_{\rm ref})$')
ax2.tick_params(axis='both')
ax2.set_xlim((1.0e-4, max(kplot)))
ax2.set_ylim((0.997, 1.059))
ax2.legend()
ax2.grid(False, which='both', axis='both')

fig.text(0.555, 0.01, r'$k ~[{\rm Mpc}^{-1}]$',
    ha='center', rotation='horizontal')

plt.savefig(rfpath+"plots/Figure_9.png") 


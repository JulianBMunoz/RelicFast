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
    "figure.dpi" : 300, 
    "figure.figsize" : [30, 30],
    'figure.subplot.left': 0.13,
    'figure.subplot.right': 0.98,
    'figure.subplot.top': 0.98,
    'figure.subplot.bottom': 0.06,
    #"figure.constrained_layout.use" : True, 
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

m_ax = np.array([
    np.power(10., -26.), 
    np.power(10., -29.), 
    np.power(10., -32.),
####################################
    np.power(10., -26.), 
    np.power(10., -29.), 
    np.power(10., -32.)
])
omega_ax = np.append(
    np.array(int(len(m_ax)/2)*[0.05*omega_cdm_LCDM]),
    np.array(int(len(m_ax)/2)*[1.0e-9*omega_cdm_LCDM])
)
sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-4
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
axion_aosc = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
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
            #"hubble = 0.701", "hubble = "+f'{h[ax_idx]:.6f}'
            "hubble = 0.701", "hubble = "+f'{h:.6f}'
        ).replace(
            "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax[ax_idx]:.6e}'
        ).replace(
            "N_klong = 1", "N_klong = "+str(Nk)
        ).replace(
            "omegac = 0.11271", "omegac = "+f'{omega_cdm[ax_idx]:.6e}'
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
    
    os.system('./relicfast run.ini')
    
    # Collect axion background evolution
    axion_background.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_background.dat"))
    axion_aosc.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_aosc.dat"))
    
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
    
    PM = []
    TF = []
    for z_idx, z_val in enumerate(z_vals): 
        PM.append(np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(len(z_vals)-z_idx)+'.dat'))
        TF.append(np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_val:.3f}'))
    
    data_pm.append(PM)
    data_tf.append(TF)
    
    #data_eulbias.append(
    #    np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    #)
    #data_lagbias.append(
    #    np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
    #)
    #print(data_lagbias[-1][0])

    os.system('mv ./run.ini ./run_'+str(ax_idx)+'.ini')
    os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(ax_idx)+'.ini')  

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

#Plot Pmm(k|z)

redshifts = np.array([0., 25., 200.])
#redshifts = np.array(z_vals)


############################
shift = int(len(m_ax)/2)

tflookup = {
    "k/h" : 0,
    "cdm" : 1,
    "baryon" : 2,
    "photon" : 3, 
    "radiation" : 4,
    "massive nu" : 5,
    "axion" : 6
}



fig, ax = plt.subplots(len(redshifts), 1,
    sharex=True,
    figsize=(30, 20.*len(redshifts)),
    gridspec_kw={'height_ratios': [1]*len(redshifts)}
)
fig.subplots_adjust(hspace=0)    

for z_idx, z_val in enumerate(redshifts): 
    pm_idx = np.argmin(np.abs(z_vals - z_val))
    print('Requested/found redshift: ', z_val, z_vals[pm_idx])
    
    colors=sns.color_palette("flare", int(len(m_ax)/2))
    plot_x = np.logspace(-4, 0, 101)

    for ax_idx, ax_val in enumerate(m_ax[0:shift]):   
        
        kfs = np.pi * np.sqrt(m_ax[ax_idx]*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+z_val, 3.), 0.5)

        #PM = data_pm[ax_idx][len(z_vals)-pm_idx-1] #In units of (Mpc/h)^3
        #PM_lcdm = data_pm[ax_idx+shift][len(z_vals)-pm_idx-1] #In units of (Mpc/h)^3
            
        #pm_interp = scipy.interpolate.interp1d(PM[:,0], PM[:,1]) # P[(Mpc/h)^3](k[h/Mpc])
        #pm_lcdm_interp = scipy.interpolate.interp1d(PM_lcdm[:,0], PM_lcdm[:,1])
        
        TF = data_tf[ax_idx][pm_idx]
        TF_lcdm = data_tf[ax_idx+shift][pm_idx]
        
        TF_k = TF[:, tflookup['k/h']]*h # Careful of units
        TF_cdm = TF[:, tflookup['cdm']]*np.power(TF_k, 2.)
        TF_b = TF[:, tflookup['baryon']]*np.power(TF_k, 2.)
        TF_x = TF[:, tflookup['axion']]*np.power(TF_k, 2.)
        
        TF_lcdm_k = TF_lcdm[:, tflookup['k/h']]*h # Careful of units 
        TF_lcdm_cdm = TF_lcdm[:, tflookup['cdm']]*np.power(TF_lcdm_k, 2.)
        TF_lcdm_b = TF_lcdm[:, tflookup['baryon']]*np.power(TF_lcdm_k, 2.)
        TF_lcdm_x = TF_lcdm[:, tflookup['axion']]*np.power(TF_lcdm_k, 2.)
        
        TF_cdm_interp = scipy.interpolate.interp1d(TF_k, TF_cdm)
        TF_b_interp = scipy.interpolate.interp1d(TF_k, TF_b)
        TF_x_interp = scipy.interpolate.interp1d(TF_k, TF_x)
        
        TF_lcdm_cdm_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_cdm)
        TF_lcdm_b_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_b)
        TF_lcdm_x_interp = scipy.interpolate.interp1d(TF_lcdm_k, TF_lcdm_x)
    
        #omega_ax[ax_idx]
        omega_M = omega_cdm[ax_idx] + omega_b_LCDM
        
        Pmm_vals = (
            Pprim(plot_x)
            * (
                np.power(f_b[ax_idx], 2.)*np.power(TF_b_interp(plot_x), 2.)
                + 2.*f_b[ax_idx]*f_cdm[ax_idx]*TF_b_interp(plot_x)*TF_cdm_interp(plot_x)
                + np.power(f_cdm[ax_idx], 2.)*np.power(TF_cdm_interp(plot_x), 2.)
            )
        )
        
        Pmm_lcdm_vals = (
            Pprim(plot_x)
            * (
                np.power(f_b[ax_idx+shift], 2.)*np.power(TF_lcdm_b_interp(plot_x), 2.)
                + 2.*f_b[ax_idx+shift]*f_cdm[ax_idx+shift]*TF_lcdm_b_interp(plot_x)*TF_lcdm_cdm_interp(plot_x)
                + np.power(f_cdm[ax_idx+shift], 2.)*np.power(TF_lcdm_cdm_interp(plot_x), 2.)
            )
        )
        
        #reco_Pmm_interp = scipy.interpolate.interp1d(TF_k, Pmm_vals*np.power(TF_k, 4.))
        #reco_Pmm_lcdm_interp = scipy.interpolate.interp1d(TF_lcdm_k, Pmm_lcdm_vals*np.power(TF_lcdm_k, 4.))
    
        plot_y =  Pmm_vals/Pmm_lcdm_vals
        #plot_y = TF_b_interp(plot_x)/TF_lcdm_b_interp(plot_x)
        #plot_y = TF_cdm_interp(plot_x)/TF_lcdm_cdm_interp(plot_x)
        
        plot_y_interp = scipy.interpolate.interp1d(np.log10(plot_x), plot_y) 

        z_osc = (1./axion_aosc[ax_idx])-1.
    
        ax[z_idx].plot(
            plot_x, 
            plot_y, 
            label=r"$m_\phi = 10^{"+f'{np.log10(ax_val):.0f}'+r"} {\rm ~eV}$",#+", $z_{osc}=$"+f"{z_osc:.2f}", 
            color=colors[ax_idx], 
            linewidth=5.
        )
        try: 
            ax[z_idx].scatter(
                kfs, 
                plot_y_interp(np.log10(kfs)), 
                facecolors='none',
                edgecolors=colors[ax_idx],
                marker='o', 
                linewidth=3.,
                s=400, 
                zorder=3
            )
        except:
            pass
        
    #plt.plot(
    #    [plot_x[0], plot_x[0]], 
    #    [plot_y[0], plot_y[0]], 
    #    color='black', 
    #    linestyle='dotted', 
    #    label="$k_{fs}$" )
    #    plt.xlabel('k')
    #    #plt.ylabel('$P_{mm}/P_{mm, \Lambda CDM}$')
    #    #plt.ylabel('$T_{b}/T_{b, \Lambda CDM}$')
    #    plt.ylabel('$T_{cdm}/T_{cdm, \Lambda CDM}$'
    #)

    ax[z_idx].set_xscale('log')
    ax[z_idx].tick_params(axis='both')
    ax[z_idx].grid(False)
    ax[z_idx].set_xlim(1.0e-4, 1.0e0) 
    ax[z_idx].set_ylim((0.65, 1.1))
    ax[z_idx].text(0.1, 1.05, r"$z = "+f"{z_val:.0f}"+r"$", 
        bbox=dict(facecolor='white', edgecolor='black', pad=10.0))

#    if z_idx==0:    
#        pass 
#    elif z_idx==1: 
#        ax[z_idx].set_ylabel(r'$P_{\rm m}/P_{\mathrm{m, }\Lambda \mathrm{CDM}}$')
#    elif z_idx==2: 
#        ax[z_idx].legend(loc='lower left')
#        ax[z_idx].set_xlabel(r"$k {\rm ~[Mpc}^{-1}{\rm ]}$") 

fig.text(0.001, 0.5, r'$P_{\rm m}/P_{\mathrm{m, }\Lambda \mathrm{CDM}}$',
    rotation='vertical', va='center')
fig.text(0.555, 0.01, r"$k {\rm ~[Mpc}^{-1}{\rm ]}$",
    ha='center', rotation='horizontal')

plt.savefig(rfpath+"plots/Figure_6.png")
    

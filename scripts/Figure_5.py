import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()
sns.set_style(style='white')

method = 3

if method==1: 
    relicfast_data_eulbias = []
    relicfast_data_lagbias = []
    axioncamb_data_tf = []
    axioncamb_data_pk = []
    relicfast_data_tf = []
    relicfast_data_pk = []
    
    rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
    rfpath_outputsuffix = "output/result-0/"
    rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"
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
    
    
    M_ax = [1.0e-30, 1.0e-26]
    
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
    
    # Calculate Bias using axionCAMB
    for m_idx, m_val in enumerate(M_ax):
        for o_idx, o_val in enumerate(omega_ax):
            print("Running RelicFast + axionCAMB for m_ax = ", M_ax[m_idx])
            
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
            
            os.system('./relicfast run.ini')
            
            relicfast_data_eulbias.append(
                np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{redshift:.2f}'+'_M'
                           +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                           +'_Nk100.dat', skiprows=1)
            )
            relicfast_data_lagbias.append(
                np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{redshift:.2f}'+'_M'
                           +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                           +'_Nk100.dat', skiprows=1)
            )
            relicfast_data_tf.append(
                np.loadtxt(rfpath_outputsuffix+'z'+f'{redshift:.2f}'+'TF_CAMB.dat', skiprows=1)
            )
            relicfast_data_pk.append(
                np.loadtxt(rfpath_outputsuffix+'power_spectra_'+'z'+f'{redshift:.2f}'+'_M'
                           +f'{np.log10(M_halo)-np.log10(h_lcdm):.2f}'
                           +'_Nk100.dat', skiprows=1)
            )
    
            #############
            # Collect axioncamb power spectrum for requested redshift 
            #output_dirs = os.listdir(rfpath_boltzmannsuffix)
            #output_dirs.remove('_params.ini')
            #
            #z_vals = np.array([])
            #
            #for str_idx, str_val in enumerate(output_dirs):
            #    if (str_val[0:15]=='_transfer_out_z'): 
            #        z_vals = np.append(
            #            z_vals, 
            #            np.float(str_val.split('_transfer_out_z')[1])
            #        )
            #            
            #z_vals = np.sort(z_vals)
            #    
        
            #pm_idx = np.argmin(np.abs(z_vals - redshift))
            #print('Requested/found redshift: ', redshift, z_vals[pm_idx])
            #
            #print('Loading: ')
            #
            #print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            #print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
            #
            #axioncamb_data_pk.append(
            #    np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
            #)
            #
            #axioncamb_data_tf.append(
            #    np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
            #)
    
            #############
            
            os.system('mv ./run.ini ./run_axioncamb_'+str(m_idx)+'.ini')
    
    
    kplot = np.geomspace(10**-3, 10**0.3, 331) #Units: h Mpc^-1
    
    colors = sns.color_palette("magma", len(omega_ax))
    
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    linestyles=["solid", "dashed"]
    for m_idx, m_val in enumerate(M_ax):
        #ref = np.loadtxt(rfpath+"scripts/2104.07802_FIG2_REF_"+str(m_idx)+".csv", delimiter=",")
        
        for o_idx, o_val in enumerate(omega_ax): 
            idx = len(omega_ax)*m_idx + o_idx
            
            if o_idx==0: 
                k_lcdm_vals = relicfast_data_pk[idx][:, 0]/h_lcdm
                pm_lcdm_vals = relicfast_data_pk[idx][:, 1]
                pm_lcdm_interp = scipy.interpolate.interp1d(np.log10(k_lcdm_vals), pm_lcdm_vals)
                pm_lcdm_plot = pm_lcdm_interp(np.log10(kplot))
            else:        
                kvals = relicfast_data_pk[idx][:, 0]/h_lcdm # h/Mpc
                pmvals = relicfast_data_pk[idx][:, 1] # Pmm
        
                pminterp = scipy.interpolate.interp1d(np.log10(kvals), pmvals)
                pmplot = pminterp(np.log10(kplot))
            
                ax.plot(
                    kplot, 
                    pmplot/pm_lcdm_plot, 
                    label=r'$\Omega_{\phi}/\Omega_{d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}', 
                    color=colors[o_idx], 
                    linestyle=linestyles[m_idx]
                )
            
            #if (o_idx==(len(omega_ax)-1)):
            #    ax[m_idx].scatter(ref[:,0][::10], ref[:,o_idx+1][::10], marker='.', color="black", label="2104.07802")
            #else: 
            #    ax[m_idx].scatter(ref[:,0][::10], ref[:,o_idx+1][::10], marker='.', color="black")
        
        ax.set_xscale('log')
        ax.set_xlabel(r'$k ~[h ~{\rm Mpc}^{-1}]$', fontsize=30)
        ax.set_ylabel(r'$P_m/P_{m, ~LCDM}$', fontsize=30)
        ax.grid(False, which='both', axis='both')
        ax.tick_params(axis='both', labelsize=25)
    
        if (m_idx==0): 
            ax.legend(fontsize=25)
    #        ax[m_idx].title((r'$M_\chi = 10^{'+f'{np.log10(m_val):.1f}'+r'}$ eV'), fontsize=30)
    
    plt.savefig("Figure_4.png")    

elif method==2: 
    axioncamb_data_pk = []

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

    Omega_ax = np.array([1e-9, 0.01, 0.04, 0.07])*omega_cdm_LCDM/np.power(h_lcdm, 2.)
    omega_ax = Omega_ax*np.power(h_lcdm, 2.)
    omega_cdm = omega_cdm_LCDM - omega_ax

    M_ax = [1.0e-30, 1.0e-26]

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

    kplot = np.geomspace(10**-3, 10**0.3, 100) #Units: [h Mpc^-1]
    colors = sns.color_palette("magma", len(omega_ax))
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    linestyles=["solid", "dashed"]
    for m_idx, m_val in enumerate(M_ax):
        #ref = np.loadtxt(rfpath+"scripts/2104.07802_FIG2_REF_"+str(m_idx)+".csv", delimiter=",")

        for o_idx, o_val in enumerate(omega_ax):
            idx = len(omega_ax)*m_idx + o_idx

            if o_idx==0:
                k_lcdm_vals = axioncamb_data_pk[idx][:, 0] #Units: [h Mpc^-1]
                pm_lcdm_vals = axioncamb_data_pk[idx][:, 1]*np.power(h_lcdm, -3.) #Units: [Mpc^-3]
                pm_lcdm_interp = scipy.interpolate.interp1d(np.log10(k_lcdm_vals), pm_lcdm_vals)
                pm_lcdm_plot = pm_lcdm_interp(np.log10(kplot))
            else:
                kvals = axioncamb_data_pk[idx][:, 0] #Units: [h Mpc^-1]
                pmvals = axioncamb_data_pk[idx][:, 1]*np.power(h_lcdm, -3.) #Units: [Mpc^-3]

                pminterp = scipy.interpolate.interp1d(np.log10(kvals), pmvals)
                pmplot = pminterp(np.log10(kplot))

                ax.plot(
                    kplot,
                    pmplot/pm_lcdm_plot,
                    label=r'$\Omega_{\phi}/\Omega_{d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}',
                    color=colors[o_idx],
                    linestyle=linestyles[m_idx]
                )

        ax.set_xscale('log')
        ax.set_xlabel(r'$k ~[h ~{\rm Mpc}^{-1}]$', fontsize=30)
        ax.set_ylabel(r'$P_m/P_{m, ~LCDM}$', fontsize=30)
        ax.grid(False, which='both', axis='both')
        ax.tick_params(axis='both', labelsize=25)

        if (m_idx==0):
            ax.legend(fontsize=25)
    #        ax[m_idx].title((r'$M_\chi = 10^{'+f'{np.log10(m_val):.1f}'+r'}$ eV'), fontsize=30)

    plt.savefig("Figure_4.png")

elif method==3: 
    
    axioncamb_data_pk = []
    axioncamb_data_tf = []

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
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
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
                    label=r'$\Omega_{\phi}/\Omega_{d}$ = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}',
                    color=colors[o_idx],
                    linestyle=linestyles[m_idx], 
                    linewidth=5.0
                )

        ax.set_xscale('log')
        ax.set_xlabel(r'$k ~[h ~{\rm Mpc}^{-1}]$', fontsize=30)
        ax.set_ylabel(r'$P_m/P_{m, ~LCDM}$', fontsize=30)
        ax.grid(False, which='both', axis='both')
        ax.tick_params(axis='both', labelsize=25)

        if (m_idx==0):
            ax.legend(fontsize=25)
    #        ax[m_idx].title((r'$M_\chi = 10^{'+f'{np.log10(m_val):.1f}'+r'}$ eV'), fontsize=30)

    plt.savefig("Figure_4.png")

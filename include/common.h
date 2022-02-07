    
    #ifndef __COMMON_H_INCLUDED__
    #define __COMMON_H_INCLUDED__
    
    
    #define _CLASS_ 0
    #define _CAMB_ 1
    #define _AXIONCAMB_ 2
    
    #define _FALSE_ 0
    #define _TRUE_ 1
    
    
    //Boltzmann code options:
    #define boltzmann_tag  _AXIONCAMB_    //0 for CLASS, 1 for CAMB. CLASS needed
    //in most cases.
    
    #define run_boltzmann_option   _TRUE_    //0 if no, otherwise YES. If
    //choosing more than one
    //z_collapse have to set to 1.
    
    //collapse code options:
    #define do_clustering_tag   _FALSE_    //whether to do clustering of
    //extra species or not.
    
    #define debug_mode  0    //whether to print out and save to file some
    //tests. 0 (or _FALSE_) for none, 1 (or _TRUE_)
    //for some, 2 or larger for A LOT.
    
    //and precision parameters:
    
    #define precision_scale 10    //from 1 to 10, precision in the
    //scale-dependence of b, through \delta_crit.
    //Recommended 2, above 4 barely any
    //difference.
    
    #define precision_normalization 4    //from 1 to 10, precision in the
    //normalization of the bias, it also
    //makes the collapse ODE solver use
    //more points, taking longer.
    //Recommended 5, above 10 barely any
    //difference.
    
    #define boost_initial_conditions 1.0    //To make initial conditions of
    //\delta_short bigger by that
    //amount. Code will complain if it
    //cannot find \delta_short.
    
    #define delta_short_constant_centroid 0.015    //the centroid of the
    //delta_shorts over which we
    //search.
    
    #define delta_long_max (delta_short_constant_centroid/10.0)    //the maximum
    //delta_long
    //we take to
    //calculate
    //the bias.
    //If you
    //change z_i
    //check that
    //it still
    //works.
    
    #define N_delta_long 2    //number of delta_longs for derivative. 2 is
    //sufficient, since it's extremely linear.
    
    #define precision_clustering 1.0    //from 0.5 to 2, makes the BKT
    //clustering integral use more points.
    //1 yields % level precision for 1-eV
    //neutrinos.
    
    #define Nrecusions_nu_clustering 1*(do_clustering_tag>0)    //how many
    //recursions of
    //finding Mnu
    //from R and
    //viceversa we
    //do.
    
    #define Mhalo_min_clustering   1e14    //in Msuns
    
    #define mnu_min_clustering 0.3    //minimum neutrino mass to consider doing
    //clustering, at M_halo =
    //Mhalo_min_clustering Msun, we scale it
    //as 1/sqrt(Mhalo); [1310.6459] multiplied
    // by T_extra/Temp_nu for extra species
    //(since it's m/T that matters). Feel free
    //to change these parameters to avoid
    //doing the clustering integral. If you
    //only care about the scale dependence of
    //b(k) you should be fine. With these
    //parameters you should still preserve
    //sub-percent-level precision in the
    //(total) bias. The effect is mostly scale
    //independent.
    
    
    //structure with all cosmological and input parameters:
    //the inputs are read with read_input_file(), and Cosmology is initialized
    //with prepare_cosmology(),
    struct Cosmology{
    
    //parameters read as input
    int tag_thermal_relic;    //whether we include thermal relics
    int tag_sterile_nu;    //whether we include sterile neutrinos
    
    
    double m_TR;    //mass of thermal relic
    double omega_TR;    //energy density of TR
    double m_SN;    //mass of sterile nu (T is assumed to be T0_nu)
    
    double omegab;    //baryon physical density (Omegab*h^2)
    double omegac;    //CDM physical density (Omegac*h^2)
    double h;    //hubble param. H0/100
    double As;    //amplitude of perturbations
    double ns;    //scalar tilt
    double mnu1;    //mass of neutrino 1
    double mnu2;    //mass of neutrino 2
    double Neff_input;    //N_effective read, we will subtract 1 per each
    //massive neutrino to get Neff.
    
    double omega_ax;    //Axion abundance. Small omega. omega_ax = h^2 Omega_ax
    double m_ax;    //Axion mass in units of eV.
    
    //derived parameters:
    
    double Neff;    //From Neff_input subtracting neutrinos that are massive
    //and so on.
    
    double omega_SN;    //sterile neutrino energy density. m_SN/94.07
    
    double omega_extra;    //energy density in TR or SN
    double Omega_extra;    //omega_extra/h^2
    double m_extra;    //mass of the extra species (TR or SN)
    
    
    int counter_massive_nus;    //how many massive nus (0, 1, or 2, since 3rd
    //nu comes in as a SN)
    
    double omeganu1;    //energy density in nu1
    double Omeganu1;
    
    double omeganu2;    //energy density in nu2
    double Omeganu2;
    
    double Omega_ax;    //energy density in axion
    
    double OmegaM;    //CMB+b energy density (NOT INCLUDING RELICS/NUs!). Used
    //for halo masses
    
    
    double T0_photons;    //temperature of photons today
    double T0_nu;    //neutrino temperature at z=0 in K. Simpler to keep
    // 1/94.07, there are more refined models.
    
    double T0_TR;    // Fixed through cosmic abundance, assuming relativistic at
    // decoupling and non-rel today
    
    double T0_extra;    // either T0_TR or TO_nu (TR or SN)
    
    
    double OmegaG;    //Omega_photons at z=0
    
    double Omeganu_massless;    //Omega_masslessneutrinos
    
    double OmegaR;    // the part that is radiation (at all z)
    
    double OmegaL;     //Omega_Lambda, we close the Friedmann Eq. assuming flat
    //Universe.
    
    double H0_Mpc;    // H0 in Mpc-1
    
    
    //collapse and other parameters:
    //read:
    
    int N_zcoll;    //how many z_collapses
    int N_Mhalo;    //how many halo masses
    int N_klong;    //how many k_longs.
    
    double z_collapse;    //the redshift of collapse for each iteration.
    double z_collapse_bot;    //minimum and maximum redshifts over which we
    //calculate. Linearly spcaed.
    
    double z_collapse_top;
    
    double Mhalo;    //the halo mass for each iteration, in Msun.
    double Mhalo_Mpc;    //G*Mhalo, in Mpc.
    double Mhalo_min;    //mhalo (in Msun) minimum and maximum over which we
    //calculate. Logspaced.
    double Mhalo_max;
    
    double ktop;    //Mpc-1, maximum and minimum k over which we calculate
    //things.
    
    double kbot;
    
    //derived
    
    int file_tag;    //tag for filenames. choose as you wish.
    
    double *z_collapse_array;    //array with z_collapses.
    int *j_collapse_array;    //which j indices correspond to which z_collapse.
    
    double *Mhalo_array;    //array with masses that we calculate for.
    
    double *klong_list_input;    //for which k_longs we calculate. We use
    //slightly different values to match CLASS/CAMB
    //for higher precision.
    
    double *axion_a;
    double *axion_z;
    double *axion_w;
    double *axion_rho;
    double *axion_p;
    double *axion_osc;
    int *axion_N;
    
    };
    
    #if boltzmann_tag==0
    #define Ninput 22    //how many inputs the code takes
    #define Ninput_int 5    //how many of the inputs are integers instead of
    //doubles.
    #elif boltzmann_tag==1
    #define Ninput 22    //how many inputs the code takes
    #define Ninput_int 5    //how many of the inputs are integers instead of
    //doubles.
    #elif boltzmann_tag==2
    #define Ninput 24    //how many inputs the code takes
    #define Ninput_int 5    //how many of the inputs are integers instead of
    //doubles.
    #endif
    
    //other cosmological constants that do not change:
    #define omega_constant  4.4806e-7    //Boltzmann factor divided by
    // rho_crit * h^2
    
    #define kpivot 0.05    //in Mpc-1, chosen to match Planck2015
    
    #define m_SN_thresh_Neff  1.0    //In eV; threshold to consider sterile
    //neutrinos as Neff or not (to subtract 1
    //from Neff).
    
    //if 3+SN have m_SN_thresh_Neff > m_SN, otherwise not.
    
    #define mass_constant_nu 93.14    //m_nu/\omega_nu in eV, good for
    //T_nu = T_\gamma * 0.71599
    
    #define ratio_T_nu_T_gamma 0.71599    //recommended by CLASS for more
    //accurately approximate neutrino
    //decoupling, if set to (4/11)^(1/3)
    //change ratio_T_nu_T_gamma to 94
    
    #define constant_change_Neff 1.0132    //this is by how much Neff changes
    //when subtracting each neutrino.
    //Not exactly one (see CLASS).
    
    //halo mass functions we implement:
    #define _MICE_ 1
    #define _wCDM_ 2
    #define _ST_ 3
    #define HMF_option _ST_    //Option for what Halo Mass Function to use.
    //1:MICE, 2:wLCDM, 3:ST. Add your own in Bias.h
    
    // Other code parameters:
    
    #define zi 200.0    //initial redshift for collapse. Same than LoVerde2014,
    //we ignore non-linear evolution until then other than
    //through CAMB.
    
    #define zf 0.0    //haloes should collapse by today--at the latest.
    
    #define length_transfer_camb 597    //Number of elements in the transfer
    //function output from CAMB
    
    #define length_transfer_axioncamb 597    //Number of elements in the transfer
    //function output from axionCAMB
    
    #define length_transfer_class 122    //Number of elements in the transfer
    //function output from CLASS
    
    #define length_transfer ( \
    (boltzmann_tag==_CLASS_)*length_transfer_class \
    + (boltzmann_tag==_CAMB_)*length_transfer_camb \
    + (boltzmann_tag==_AXIONCAMB_)*length_transfer_axioncamb \
    )
    
    #define Nz_transfer 100    //how many redshifts we take for transfer
    //functions. For more than ~115 need to skip
    //lines in CLASS file, so it's not supported.
    
    #define print_headers _TRUE_    //whether to print headers in output files
    //or not
    
    #define z_collapse_min 0.001    //smallest z_collapse before the code
    //complains that it has to change the
    //binning to linear instead of logarithmic.
    
    ///Mathematical constants and unit conversions:
    
    #define PI 3.1415926535
    
    #define Msuntokm 1.47578    //Msun to km
    #define kmtoMpc  3.24e-20    //km to Mpc
    #define MsuntoMpc 4.80e-20    //Msun to Mpc.
    #define c_light  299792    //speed of light in km/s
    #define pctoly  3.262    //parsec to lightyear
    
    #define KtoeV 8.617e-5    //we take k_B=1, so we convert K to eV and
    //viceversa.
    
    #define  hbareVs 6.582e-16    // hbar in eV * s, we will assume hbar=1
    //throughout, so we have to use to convert
    //back sometimes.
    
    #define  hbareVkm  1.975e-10    // hbar in eV * km
    #define  hbareVMpc  6.398e-30    // hbar in eV * Mpc
    
    #define  eVtokg    1.783e-36    // 1 eV in kg
    #define  Msuntokg   1.9886e30    // 1 Msun in kg
    #define  eVtoMpc   4.304e-86    //eVtokg/Msuntokg*MsuntoMpc; here we assume
    //G=1: geometrized units.
    
    #endif

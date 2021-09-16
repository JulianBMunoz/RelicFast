////////////////////////////////////////////////////////////////////////////////////////
//
//    Code to find the bias of a halo of mass Mhalo
//    Input \delta_initial and \delta_crit
//        output b_L and b_E (Lagrangian and Eulerian biases)
//    By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////


#include "bias.h"

int get_bias(Cosmology *cosmo, double *zlist_transfer){

    int lengthname=200;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;
    double *klong_list;//in 1/Mpc. Input and closest CAMB values.
    klong_list = allocate_1D_array(cosmo->N_klong);
    int i_delta_long, i_klong, i, j;
    double **delta_short_crit; //delta_short_critical, 
                               //extrapolated to z_collapse.
    double **delta_long_collapse; //delta_long, extrapolated to z_collapse.

    delta_short_crit=allocate_2D_array(N_delta_long,cosmo->N_klong);
    delta_long_collapse=allocate_2D_array(N_delta_long,cosmo->N_klong);

    char *buffer; //To save header and throw out
    buffer=(char *)malloc((lengthname+1)*sizeof(char));

    //which z_collapse are we in:
    int number_collapse=find_value(
        cosmo->N_zcoll,
        cosmo->z_collapse_array,
        cosmo->z_collapse
    );
    int j_collapse = cosmo->j_collapse_array[number_collapse];
    if (debug_mode>0) printf(
        "(in bias) number_collapse=%d, j_collapse=%d \n",
        number_collapse,
        j_collapse
    );

    //we will loop over Mhalo and klong, to avoid reading transfer functions 
    //many time. there isn't a lot to win from looping over z inside here, as 
    //opposed to outside, since you still have to read Boltzmann output every 
    //time.

    /////////////////////////////////////////
    //// we read the transfer functions    
    ////////////////////////////////////////

    double *k_transfer_array, *transfer_array_z_coll; //for CDM+b, derivative
    double *transfer_array_gamma_z_coll;
    double *transfer_array_nu_massless_z_coll;
    double *transfer_array_nu1_z_coll; 
    //for photons, massless neutrinos, and massive neutrino1

    k_transfer_array = allocate_1D_array(length_transfer);
    transfer_array_z_coll = allocate_1D_array(length_transfer);
    transfer_array_gamma_z_coll = allocate_1D_array(length_transfer);
    transfer_array_nu_massless_z_coll = allocate_1D_array(length_transfer);
    transfer_array_nu1_z_coll = allocate_1D_array(length_transfer);

    double *transfer_array_nu2_z_coll, *transfer_array_extra_z_coll; 
    //for mnu2 and extra species

    transfer_array_nu2_z_coll = allocate_1D_array(length_transfer);
    transfer_array_extra_z_coll = allocate_1D_array(length_transfer);

    int transfer_check;

    if(boltzmann_tag == 0){//CLASS
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/z%d_tk.dat",
            boltzmann_tag,
            cosmo->file_tag, 
            j_collapse+1
        ); 
        //zlist reversed, starts at 1
    }
    else { //CAMB
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",
            boltzmann_tag, 
            cosmo->file_tag, 
            cosmo->z_collapse
        ); 
        //We reuse the same filename variable name.
    }
    transfer_check=gettransfer_matter(
        cosmo, 
        filename, 
        k_transfer_array, 
        transfer_array_z_coll
    );

    do_check(transfer_check);

    //we have a different function to extract photon transfer function
    transfer_check=gettransfer_gamma(cosmo, filename, k_transfer_array, transfer_array_gamma_z_coll);
    do_check(transfer_check);
    //same for massless neutrinos
    transfer_check=gettransfer_nu_massless(
        cosmo, 
        filename, 
        k_transfer_array, 
        transfer_array_nu_massless_z_coll
    );
    do_check(transfer_check);
    //and for massive neutrino 1
    if(cosmo->mnu1>0){
        transfer_check=gettransfer_nu1(
            cosmo, 
            filename, 
            k_transfer_array, 
            transfer_array_nu1_z_coll
        );
        do_check(transfer_check);
    }
    //and for massive neutrino 2
    if(cosmo->mnu2>0){
        transfer_check=gettransfer_nu2(
            cosmo, 
            filename, 
            k_transfer_array, 
            transfer_array_nu2_z_coll
        );
        do_check(transfer_check);
    }
    //and for extra species
    if(cosmo->Omega_extra > 0){
        transfer_check=gettransfer_extra(
            cosmo, 
            filename, 
            k_transfer_array, 
            transfer_array_extra_z_coll
        );
        do_check(transfer_check);
    }

    //we also save the transfer function from CLASS/CAMB at z_collapse in the 
    //output folder.
    if(boltzmann_tag == 0){//CLASS
        lengthname=sprintf(
            filename,
            "cp Boltzmann_%d/transfer_files_%d/z%d_tk.dat \
                output/result-%d/z%.2ftk.dat",
            boltzmann_tag, 
            cosmo->file_tag, 
            j_collapse+1, 
            cosmo->file_tag, 
            cosmo->z_collapse
        ); //zlist reversed, starts at 1
    }
    else { //CAMB
        lengthname=sprintf(
            filename,
            "cp Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f \
                output/result-%d/z%.2fTF_CAMB.dat",
            boltzmann_tag, 
            cosmo->file_tag, 
            cosmo->z_collapse, 
            cosmo->file_tag, 
            cosmo->z_collapse
        ); 
        //We reuse the same filename variable name.
    }
    system(filename);

    //and we will compute sigma_M(z_coll).
    double sigmaM_coll;

    /////////////////////////////////////////
    //// these are the HMFs that we use    ////
    ////////////////////////////////////////

    double HMF;
    //choose your favorite halo mass function:
    //the input (called HMF) is simply 
    //dlogn/d delta_crit = d log f/d delta_crit.
    //if your HMF does not have an analytic derivative do numerical instead.

    //Option 1: MICE (0907.0019)
    //the fit is for f_MICE (\sigma), so we transform 
    //\sigma -> \sigma/\delta_crit and take the derivative.

    // HMF = -2 c delta_crit/(delta_ref^2 * sigma^2) 
    // + a/delta_crit/(1+b*(delta_ref * sigma/delta_crit)^a), 
    //where delta_ref=1.686, since they only quote it as a function of \sigma.


    double delta_ref=1.686;
    double a_MICE=1.37*pow(1+cosmo->z_collapse,-0.15);
    double b_MICE=0.3*pow(1+cosmo->z_collapse,-0.084);
    double c_MICE=1.036*pow(1+cosmo->z_collapse,-0.024);

    double HMF_MICE;

    //Option 2: arXiv:1005.2239
    // HMF =d log(n)/d \delta_crit = (q-a (delta_crit/sigma)^2)/delta_crit 
    // - (2 p/delta_crit)/(1+(a delta_crit^2/sigma^2)^p)

    double q_hmf=1.795;
    double p_hmf=0.807;
    double a_hmf_0=0.788; // a_hmf (z) = a_hmf_0/(1+z)^a_hmf_exp
    double a_hmf_exp=0.01;
    double a_hmf=a_hmf_0/pow(1.+cosmo->z_collapse,a_hmf_exp);

    double HMF_Batt;

    //Option 3: Sheth Tormen (astro-ph/9901122)
    // HMF = (1-a (delta_crit/sigma)^2)/delta_crit 
    // - (2 p/delta_crit)/(1+(a delta_crit^2/sigma^2)^p)
    // it's essentially the same as Option 2, but with no z-dependence and 
    // with q set to unity;

    double a_ST=0.707; 
    //updated from original 0.707, from Ref.1005.2239 Table 3

    double p_ST=0.3;
    double HMF_ST;

    double delta_crit_sim=1.686; 
    //usual value of delta_crit used in simulations to obtain the HMF.

    double k_long; 
    //k_long for each iteration

    //lagrangian and Eulerian biases:
    //b_L = d log(n)/d \delta_crit * d \delta_crit / d \delta_long
    //b_E = Phm/Pmm.

    double *b_L, *b_E;
    b_L = allocate_1D_array(cosmo->N_klong);
    b_E = allocate_1D_array(cosmo->N_klong);

    double derivative;
    int counter; //to get average of delta_crit.

    const int N_species = (
        1 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0)
    );
    //how many species (cdm+b) + others

    int species_counter=0;

    const double Omegatot = (
        cosmo->OmegaM + cosmo->Omeganu1 + cosmo->Omeganu2 + cosmo->Omega_extra
    ); 
    //remember we call OmegaM=Omegab+Omegac

    double *fraction_species; //fractions of matter today in each species.
    fraction_species=allocate_1D_array(N_species); // cdm, nu1, nu2, extra.

    fraction_species[0] = cosmo->OmegaM/Omegatot;//CDM+b
    if(cosmo->mnu1>0){
        species_counter++;
        fraction_species[species_counter] = cosmo->Omeganu1/Omegatot;
    }
    if(cosmo->mnu2>0){
        species_counter++;
        fraction_species[species_counter] = cosmo->Omeganu2/Omegatot;
    }
    if(cosmo->Omega_extra>0){
        species_counter++;
        fraction_species[species_counter] = cosmo->Omega_extra/Omegatot;
    }
    for(species_counter=0;species_counter<N_species;species_counter++){
        if (debug_mode>0) printf(
            "fraction_%d=%.3le \n\n\n",
            species_counter,
            fraction_species[species_counter]
        );
    }
    species_counter=0; //we reset the counter.

    double *tf_species;
    tf_species=allocate_1D_array(N_species); // cdm, nu1, nu2, extra.

    double **Power_spectra; //P_ij, for i and j species. At fixed k_long

    Power_spectra=allocate_2D_array(N_species,N_species); 
    // cdm, nu1, nu2, extra.

    double *Pmm, *Pmh, *Phh; //Power spectra at different k_longs
    Pmm=allocate_1D_array(cosmo->N_klong);
    Pmh=allocate_1D_array(cosmo->N_klong);
    Phh=allocate_1D_array(cosmo->N_klong);

    double factor; //(2pi^2) and other factors.

    int iM;

    for(iM=0;iM<cosmo->N_Mhalo;iM++){
        cosmo->Mhalo = cosmo->Mhalo_array[iM];
        cosmo->Mhalo_Mpc = cosmo->Mhalo * MsuntoMpc;
        //Mhalo in Mpc (*G)

        /////////////////////////////////////////////
        //// we read the results from files.    ////
        ///////////////////////////////////////////

        lengthname=sprintf(
            filename,
            "output/result-%d/delta_crit_z%.2f_M%.2f_Nk%d.dat",
            cosmo->file_tag,cosmo->z_collapse, 
            log10(cosmo->Mhalo),
            cosmo->N_klong
        ); 
        //we save the delta_crit extrapolated to z=0.5

        fp=fopen(filename,"r");

        if(print_headers!=0){
            fgets(buffer, 100, fp); 
            //reads 100 characters in first line. To get rid of headers.
        }
        for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
            for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
                fscanf(
                    fp,
                    "%le %le %le", 
                    &delta_long_collapse[i_delta_long][i_klong], 
                    &klong_list[i_klong], 
                    &delta_short_crit[i_delta_long][i_klong]
                );
            }
        }

        fclose(fp);

        sigmaM_coll=getsigma_M(cosmo, k_transfer_array, transfer_array_z_coll);

        /////////////////////////////////////
        //// we now compute the bias.    ////
        ////////////////////////////////////

        //we find the "average" value for delta_crit_sim for our collapse. 
        //Just to get a fair estimate of the amplitude of b_1^L

        delta_crit_sim=0.;
        counter=0;
        for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
            for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
                counter++;
                delta_crit_sim+=(
                    delta_long_collapse[i_delta_long][i_klong]
                    + delta_short_crit[i_delta_long][i_klong]
                );
            }
        }

        delta_crit_sim/=(1.*counter);
        if(debug_mode>0){
            printf("delta_sim=1.686, delta_us=%.3le \n",delta_crit_sim);
        }

        HMF_MICE = (
            -2. 
            * c_MICE 
            * delta_crit_sim
            / (delta_ref*delta_ref * sigmaM_coll * sigmaM_coll) 
            + a_MICE/delta_crit_sim
            / (1. + b_MICE*pow(sigmaM_coll*delta_ref/delta_crit_sim,a_MICE))
        );

        HMF_Batt= (
            (q_hmf-a_hmf*pow(delta_crit_sim/sigmaM_coll,2.))
            / delta_crit_sim 
            - (2*p_hmf/delta_crit_sim)
            / (1.+pow(a_hmf*pow(delta_crit_sim/sigmaM_coll,2.),p_hmf))
        );

        HMF_ST= (
            (1-a_ST*pow(delta_crit_sim/sigmaM_coll,2.))
            / delta_crit_sim 
            - (2*p_ST/delta_crit_sim)
            / (1.+pow(a_ST*pow(delta_crit_sim/sigmaM_coll,2.),p_ST))
        );

        //select a halo mass function, or add your own:
        if(HMF_option==1){ //MICE
            HMF=HMF_MICE;
        }
        else if(HMF_option==2){ //wLCDM
            HMF=HMF_Batt;
        }
        else if(HMF_option==3){ //Sheth-Tormen
            HMF=HMF_ST;
        }
        else {
            printf("Halo Mass Function not chosen, set HMF=HMF_MICE to 1, 2, \
                or 3. MICE selected by default \n");
            HMF=HMF_MICE;
        }

        //printf("HMF = %le \n",HMF);

        for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
            derivative=(
                (
                    delta_short_crit[N_delta_long-1][i_klong]
                    - delta_short_crit[0][i_klong]
                ) /
                (
                    delta_long_collapse[N_delta_long-1][i_klong]
                    - delta_long_collapse[0][i_klong]
                )
            ); 
            // ddelta_crit/ddelta_long

            b_L[i_klong] = HMF * derivative;
        }


        // for(i_klong=0;i_klong<length_transfer;i_klong++){
        //     printf("%le \n",k_transfer_array[i_klong]);
        // }


        for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
            k_long=klong_list[i_klong];
            factor=(
                (2.*PI*PI)
                * cosmo->As
                * pow(k_long/kpivot,cosmo->ns-1.)
                / k_long/k_long/k_long
            );
            Pmm[i_klong] = Pmh[i_klong] = Phh[i_klong] = 0.;

            tf_species[0] = interpol(
                transfer_array_z_coll,
                k_transfer_array,
                length_transfer,k_long
            );

            if(cosmo->mnu1>0){
                species_counter++;
                tf_species[species_counter] = interpol(
                    transfer_array_nu1_z_coll,
                    k_transfer_array,
                    length_transfer,
                    k_long
                );
            }
            if(cosmo->mnu2>0){
                species_counter++;
                tf_species[species_counter] = interpol(
                    transfer_array_nu2_z_coll,
                    k_transfer_array,
                    length_transfer,
                    k_long
                );
            }
            if(cosmo->Omega_extra>0){
                species_counter++;
                tf_species[species_counter] = interpol(
                    transfer_array_extra_z_coll,
                    k_transfer_array,
                    length_transfer,
                    k_long
                );
            }
            species_counter=0; //we reset the counter.

            for(i=0;i<N_species;i++){
                for(j=0;j<N_species;j++){
                    Power_spectra[i][j] = (
                        factor 
                        * tf_species[i] 
                        * tf_species[j]
                    ); // P_ij

                    Pmm[i_klong]+= (
                        Power_spectra[i][j] 
                        * fraction_species[i] 
                        * fraction_species[j]
                    ); //P_mm
                    //printf("k/h=%.3le, i=%d, j=%d, P=%.1le \n",k_long/h,i,j,Power_spectra[i][j]);
                }
                Pmh[i_klong]+= (
                    (1.0+b_L[i_klong]) 
                    * Power_spectra[i][0] 
                    * fraction_species[i]
                ); //P_mhalo
            }
            Phh[i_klong] += (
                (1.0+b_L[i_klong]) 
                * (1.0+b_L[i_klong]) 
                * Power_spectra[0][0]
            ); //P_hh = (1+bL)^2Pcc

            b_E[i_klong]=Pmh[i_klong]/Pmm[i_klong];

            if(N_species>1 && debug_mode>0){
                printf(
                    "k/h=%.1le,  Tx/Tc=%.3le, b_E=%.3le, b_L=%.3le \n",
                    k_long/cosmo->h, 
                    Power_spectra[1][0]/Power_spectra[0][0], 
                    Pmh[i_klong]/Pmm[i_klong], 
                    b_L[i_klong]
                );
            }
        }

        /////////////////////////////////////
        //// we save the output in file  ////
        ////////////////////////////////////

        lengthname=sprintf(
            filename,
            "output/result-%d/bias_Lagrangian_z%.2f_M%.2f_Nk%d.dat",
            cosmo->file_tag,
            cosmo->z_collapse, 
            log10(cosmo->Mhalo),
            cosmo->N_klong
        ); 
        //we save the delta_crit extrapolated to z=0.5

        fp=fopen(filename,"w");

        if(print_headers!=0){
            fprintf(fp,"k_long[1/Mpc]   bias Lagrangian \n");
        }
        for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
            k_long=klong_list[i_klong];
            fprintf(fp,"%le   %le \n", k_long, b_L[i_klong]);
        }
        fclose(fp);

        lengthname=sprintf(
            filename,
            "output/result-%d/bias_Euler_z%.2f_M%.2f_Nk%d.dat", 
            cosmo->file_tag, 
            cosmo->z_collapse, 
            log10(cosmo->Mhalo),
            cosmo->N_klong
        ); 
        //we save the delta_crit extrapolated to z=0.5

        fp=fopen(filename,"w");

        if(print_headers!=0){
            fprintf(fp,"k_long[1/Mpc]   bias Eulerian \n");
        }
        for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
            k_long=klong_list[i_klong];
            fprintf(fp,"%le   %le \n", k_long, b_E[i_klong]);
        }
        fclose(fp);

        lengthname=sprintf(
            filename,
            "output/result-%d/power_spectra_z%.2f_M%.2f_Nk%d.dat", 
            cosmo->file_tag, 
            cosmo->z_collapse, 
            log10(cosmo->Mhalo),
            cosmo->N_klong
        ); 
        //we save the delta_crit extrapolated to z=0.5

        fp=fopen(filename,"w");

        if(print_headers!=0){
            fprintf(fp,"k_long[1/Mpc]\t Pmm\t Pmh\t Phh \n");
        }

        for(i_klong=0;i_klong<cosmo->N_klong;i_klong++){
            k_long=klong_list[i_klong];
            fprintf(
                fp,
                "%.4le\t %.4le\t %.4le\t %.4le \n", 
                k_long, 
                Pmm[i_klong], 
                Pmh[i_klong], 
                Phh[i_klong]
            );
        }
        fclose(fp);
    }//end of Mhalo loop


    //////////////////////////////////////////
    //// we free the allocated memory /////
    ////////////////////////////////////////


    free(b_L);
    free(b_E);

    free(klong_list);

    free_2D_array(delta_short_crit,N_delta_long);
    free_2D_array(delta_long_collapse,N_delta_long);


    free(buffer);

    free(filename);

    free(k_transfer_array);


    free(transfer_array_z_coll);
    free(transfer_array_gamma_z_coll);
    free(transfer_array_nu_massless_z_coll);
    free(transfer_array_nu1_z_coll);

    free(transfer_array_nu2_z_coll);
    free(transfer_array_extra_z_coll);


    free(Pmm);
    free(Phh);
    free(Pmh);

    free_2D_array(Power_spectra, N_species);

    free(fraction_species);
    free(tf_species);

    return 1;

}

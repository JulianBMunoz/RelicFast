////////////////////////////////////////////////////////////////////////////////////////
//
//    Code to call a Boltzmann solve (CLASS or CAMB)
//    as well as to read the transfer functions from it.
//    By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////


#include "boltzmann_Calls.h"


///////////////////////////////
///We do some tests pre run///
/////////////////////////////
int tests_pre_run(Cosmology *cosmo){
    if(boltzmann_tag == _CAMB_){//if camb:
        if(cosmo->counter_massive_nus>1){
            printf("Cannot do more than one massive neutrino with CAMB, switch to CLASS \n");
            return -1;
        }
        else if(cosmo->counter_massive_nus == 1 && cosmo->tag_sterile_nu 
            + cosmo->tag_thermal_relic != 0){
            printf("Cannot do massive neutrinos and pHDM with CAMB, switch to CLASS \n");
            return -1;
        }
    }
    if(cosmo->tag_sterile_nu + cosmo->tag_thermal_relic > 1){
        printf("Do one pHDM case at a time. \n");
        return -1;
    }

    //some coding tests:
    double mnu_low_limit=0.009;
    //eV, how light can nus be before they are (somewhat) relativistic at z=0
    //(T0~0.1meV, so for 0.1 eV we are safe).

    if(cosmo->counter_massive_nus>0){
        if(cosmo->mnu1>0 && cosmo->mnu1<mnu_low_limit){
            printf("Error in neutrino specifications, nu1 has to be non-relativistic at z=0. \n");
            return -1;
        }
        if(cosmo->mnu2>0 && cosmo->mnu2<mnu_low_limit){
            printf("Error in neutrino specifications, nu2 has to be non-relativistic at z=0. \n");
            return -1;
        }
        if(cosmo->Neff + cosmo->counter_massive_nus > 4){
            printf("Are you sure you want Neff~3 AND massive neutrinos? \
                If so change YHe=0.24 in CLASS. \n");
        }
    }
    if(cosmo->tag_sterile_nu != 0){
        if(cosmo->m_SN<mnu_low_limit){
            printf("Error in sterile neutrino specifications, \
                particle has to be non-relativistic at z=0. \n");
            return -1;
        }
    }

    double z_rel; //rough redshift at which thermal relic becomes relativistic

    if(cosmo->tag_thermal_relic != 0){
        if(cosmo->m_TR<cosmo->T0_TR*KtoeV){
            printf("Error in thermal relic specifications, particle has to be \
                non-relativistic at z=0. \n");
            return -1;
        }
        z_rel = cosmo->m_TR/(cosmo->T0_TR*KtoeV);

        if(pow(cosmo->T0_TR/cosmo->T0_nu,3.)*cosmo->m_TR
            /(cosmo->T0_nu*KtoeV)/z_rel > 1.){
            printf("Thermal relic will contribute to Neff >~ 1 at high redshift. \n");
        }
    }

    printf("Number of massive neutrinos: %d \n",cosmo->counter_massive_nus);
    printf("Additional Species? %d \n",(cosmo->Omega_extra>0));
    if(cosmo->tag_sterile_nu){
        printf("Sterile Neutrino (or 3rd neutrino), with m=%2lf eV \n",cosmo->m_SN);
    }
    if(cosmo->tag_thermal_relic){
        printf("Thermal Relic, with m=%.2lf eV and T_0=%.2lf K \n", cosmo->m_TR, cosmo->T0_TR);
    }

    if(cosmo->z_collapse_top > 1. && cosmo->N_zcoll>1 && HMF_option!=_ST_){
        printf("For z>1 HMFs from simulations might not properly calibrated, \
            change to Sheth-Tormen HMF. \n");
    }

    return 1;
}


////////////////////////////////////////////////////////////////////////////////////
//Read the transfer functions from the chosen Boltzmann code
//
//input fp(file(name)) and length_transfer of file, saves result in kgrid and TF
///////////////////////////////////////////////////////////////////////////////////////
int gettransfer_matter(Cosmology *cosmo, char *filename,  double *kgrid, double *TF){

    int j;

    double temp[length_transfer];

    double TF_grid_b[length_transfer];
    double TF_grid_c[length_transfer];
    double koverh[length_transfer];

    FILE *fp2=fopen(filename,"r");

    int Ncolum;
    int length_buffer=1000;

    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));

    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file

        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line. To get rid of headers 
        //from CLASS

        for(j=0; j<length_transfer; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            koverh[j] = transfer_temp[0];
            TF_grid_b[j] = transfer_temp[2];
            TF_grid_c[j] = transfer_temp[3];    //printf("%le \n",koverh[j]);
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that

        for(j=0;j<length_transfer;++j){
            kgrid[j] = cosmo->h * koverh[j];
            TF[j] = (
                cosmo->omegab/(cosmo->omegab+cosmo->omegac) * TF_grid_b[j] 
                + cosmo->omegac/(cosmo->omegab+cosmo->omegac) * TF_grid_c[j]
            );
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB
        Ncolum=13;    
        //number of columns in CAMB output file

        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                &koverh[j], 
                &TF_grid_c[j],
                &TF_grid_b[j],
                temp,
                temp,
                temp,
                temp, 
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            kgrid[j] = cosmo->h * koverh[j];
            TF[j] = kgrid[j]*kgrid[j]*(
                (cosmo->omegac/(cosmo->omegab+cosmo->omegac))*TF_grid_c[j]
                + (cosmo->omegab/(cosmo->omegab+cosmo->omegac))*TF_grid_b[j]
            );; 
            //for CAMB we compute total matter (c+b).
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //axionCAMB
        Ncolum=9;
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le", 
                &koverh[j], 
                &TF_grid_c[j], 
                &TF_grid_b[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        ){};

        fclose(fp2);
        for(j=0;j<length_transfer;++j){
            kgrid[j] = cosmo->h * koverh[j];
            TF[j] = kgrid[j]*kgrid[j]*(
                (cosmo->omegac/(cosmo->omegab+cosmo->omegac))*TF_grid_c[j]
                + (cosmo->omegab/(cosmo->omegab+cosmo->omegac))*TF_grid_b[j]
            );
            //printf("om_c = %le \n",(cosmo->omegac/(cosmo->omegab+cosmo->omegac)));  
        }
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be either 0 or 1 \n");
        return -1;
    }
    free(buffer);
    free(transfer_temp);
    return 1;
}

int gettransfer_gamma(Cosmology *cosmo, char *filename, double *kgrid, 
    double *TF){
    //this function extracts the transfer function of photons, make sure that 
    //the column is right if you use a different CLASS/CAMB version.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int length_buffer=1000;

    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));

    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file

        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line. To get rid of headers 
        //from CLASS

        for(j=0; j<length_transfer_class; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            TF_grid[j] = transfer_temp[1];
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = TF_grid[j];
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB
        Ncolum=13; //number of columns in CAMB output file
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //axionCAMB

        Ncolum=9;
            
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );//CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that //
        for(j=0;j<length_transfer;++j){
                TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be either 0, 1, or 2 \n");
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}

int gettransfer_nu_massless(Cosmology *cosmo, char *filename,  double *kgrid, 
    double *TF){
    //this function extracts the transfer function of massless neutrinos, make
    //sure that the column is right if you use a different CLASS/CAMB version.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int length_buffer=1000;
    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));
    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file
        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line. To get rid of headers
        //from CLASS

        for(j=0; j<length_transfer_class; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            TF_grid[j] = transfer_temp[4];
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = TF_grid[j];
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB
        Ncolum=13; //number of columns in CAMB output file
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp,
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that //

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //axionCAMB
        Ncolum=9;
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be either 0 or 1 \n");
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}


int gettransfer_nu1(Cosmology *cosmo, char *filename, double *kgrid, 
    double *TF){
    //this function extracts the transfer function of neutrino 1, make sure that
    //the column is right if you use a different CLASS/CAMB version.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int colum_extract=5;
    int length_buffer=1000;
    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));
    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file
        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line. 
        //To get rid of headers from CLASS

        for(j=0; j<length_transfer_class; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            TF_grid[j] = transfer_temp[colum_extract];
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 
        for(j=0;j<length_transfer;++j){
            TF[j] = TF_grid[j];
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB
        Ncolum=13; //number of columns in CAMB output file
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp,
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //axionCAMB
        Ncolum=9;

        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp,   
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );//CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that //
        for(j=0;j<length_transfer;++j){
                TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else{
        printf(
            "Error in gettransfer, boltzmann_tag has to be either 0, 1, or 2 \n"
        );
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}

int gettransfer_nu2(Cosmology *cosmo, char *filename, double *kgrid, 
    double *TF){
    //this function extracts the transfer function of 2nd neutrino, make sure
    //that the column is right if you use a different CLASS/CAMB version.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int length_buffer=1000;
    int colum_extract=6;

    // if(boltzmann_tag==1){ //CAMB
    //     printf("Taking CAMB and more than 1 nu, all assumed identical \n");
    // }

    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));
    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file

        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line. 
        //To get rid of headers from CLASS

        for(j=0; j<length_transfer_class; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            TF_grid[j] = transfer_temp[colum_extract];
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = TF_grid[j];
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB. WE JUST COPY NEUTRINO 1
        Ncolum=13; //number of columns in CAMB output file
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp,
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that //

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //CAMB. WE JUST COPY NEUTRINO 1
        Ncolum=9;
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le", 
                temp, 
                temp, 
                temp, 
                temp, 
                temp, 
                &TF_grid[j], 
                temp, 
                temp, 
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be \
            either 0 or 1 \n");
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}


int gettransfer_extra(Cosmology *cosmo, char *filename, double *kgrid, 
    double *TF){
    //this function extracts the transfer function of extra species, make sure
    //that the column is right if you use a different CLASS/CAMB version.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int length_buffer=1000;
    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));
    int index_extra;
    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_CLASS_){ //CLASS
        Ncolum = 8 + cosmo->counter_massive_nus + (cosmo->Omega_extra>0); 
        //number of columns in CLASS output file

        index_extra = 6 + (cosmo->mnu1>0) + (cosmo->mnu2>0) - 1; 
        //7 or 8 depending on whether nu2 are massless (nu1 acts as massless
        //neutrinos anyway)

        fgets(buffer, length_buffer, fp2); 
        //reads length_buffer characters in first line.      
        //To get rid of headers from CLASS

        for(j=0; j<length_transfer_class; j++){
            for (i=0;i<Ncolum;i++){
                fscanf(fp2, "%le" , &transfer_temp[i]);
            }
            TF_grid[j] = transfer_temp[index_extra];
        }

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not.
        //We get k_grid from the rest of files to avoid issues with Ncolum.

        for(j=0;j<length_transfer;++j){
            TF[j] = TF_grid[j];
            //in class it's not over k^2
        }
    }
    else if(boltzmann_tag==_CAMB_){ //CAMB
        Ncolum=13; //number of columns in CAMB output file
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le %le %le %le %le", 
                temp, 
                temp,
                temp, 
                temp,
                temp, 
                &TF_grid[j],
                temp,
                temp,
                temp,
                temp,
                temp,
                temp,
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }
    }
    else if(boltzmann_tag==_AXIONCAMB_){ //CAMB
        Ncolum=9; 
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le",
                temp, 
                temp, 
                temp, 
                temp, 
                temp, 
                &TF_grid[j], 
                temp, 
                temp, 
                temp
            )==Ncolum;
            ++j
        );
        //CAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }   
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be either 0 or 1 \n");
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}



int gettransfer_axion(Cosmology *cosmo, char *filename, double *kgrid, 
    double *TF){
    //this function extracts the transfer function of axion species.
    //It only works if using the axionCAMB solver.

    int j;
    double temp[length_transfer];
    double TF_grid[length_transfer];
    FILE *fp2=fopen(filename,"r");
    int Ncolum;
    int length_buffer=1000;
    char *buffer; //To save header and throw out
    buffer=(char *)malloc((length_buffer+1)*sizeof(char));
    int index_extra;
    int i;
    double *transfer_temp; //to save values of transfer
    int length_temp=20;
    transfer_temp = allocate_1D_array(length_temp); 
    //less than 20 columns always

    if(boltzmann_tag==_AXIONCAMB_){
        Ncolum=9; 
        for(
            j=0;
            fscanf(
                fp2,
                "%le %le %le %le %le %le %le %le %le",
                temp, 
                temp, 
                temp, 
                temp, 
                temp,
                temp, 
                &TF_grid[j], 
                temp, 
                temp
            )==Ncolum;
            ++j
        );
        //AXIONCAMB files have Ncolum columns

        fclose(fp2);

        //We need to give the right values to TF, evaluated over k and not k/h,
        //we define kgrid for that 

        for(j=0;j<length_transfer;++j){
            TF[j] = kgrid[j]*kgrid[j]*TF_grid[j];
        }   
    }
    else{
        printf("Error in gettransfer, boltzmann_tag has to be either 0, 1, or 2 \n");
        return -1;
    }

    free(buffer);
    free(transfer_temp);

    return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////
//        This function reads a bunch of transfer files in transfer_files/ and saves them ///
//////////////////////////////////////////////////////////////////////////////////////////

int save_transfers_camb(
    Cosmology *cosmo, 
    double zlist_transfer[], 
    double *kgrid,
    double **TFm, 
    double **TFgamma, 
    double **TFnu, 
    double **TFnu_mass, 
    double **TFnu2, 
    double **TFextra
    ){
    //we will save the TF from CAMB
    //nu2 and extra are set to be identical to nu1 if all are active.

    int lengthname=200, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;
    int checkm, checknu, checknu_mass, checkgamma, checkaxion;

    for(jz=0;jz<Nz_transfer;jz++){
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f", 
            boltzmann_tag, 
            cosmo->file_tag, 
            zlist_transfer[jz]
        );

        fp=fopen(filename, "r");
        checkm=gettransfer_matter(cosmo, filename, kgrid, TFm[jz]);
        do_check(checkm);
        checkgamma=gettransfer_gamma(cosmo, filename, kgrid, TFgamma[jz]);
        do_check(checkgamma);
        checknu=gettransfer_nu_massless(cosmo, filename,  kgrid, TFnu[jz]);
        do_check(checknu);
        checknu_mass=gettransfer_nu1(cosmo, filename, kgrid, TFnu_mass[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_nu2(cosmo, filename,  kgrid, TFnu2[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_extra(cosmo, filename, kgrid, TFextra[jz]);
        do_check(checknu_mass);
        fclose(fp);
    }

    free(filename);

    return 1;
}

int save_transfers_axioncamb(
    Cosmology *cosmo,
    double zlist_transfer[],
    double *kgrid,
    double **TFm,
    double **TFgamma,
    double **TFnu,
    double **TFnu_mass,
    double **TFnu2,
    double **TFextra,
    double **TFaxion
    ){
    //we will save the TF from axionCAMB
    //nu2 and extra are set to be identical to nu1 if all are active.

    int lengthname=200, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;
    int checkm, checknu, checknu_mass, checkgamma, checkaxion;
    
    for(jz=0;jz<Nz_transfer;jz++){
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",
            boltzmann_tag,
            cosmo->file_tag,
            zlist_transfer[jz]
        );

        fp=fopen(filename, "r");
        checkm=gettransfer_matter(cosmo, filename, kgrid, TFm[jz]);
        do_check(checkm);
        checkgamma=gettransfer_gamma(cosmo, filename, kgrid, TFgamma[jz]);
        do_check(checkgamma);
        checknu=gettransfer_nu_massless(cosmo, filename,  kgrid, TFnu[jz]);
        do_check(checknu);
        checknu_mass=gettransfer_nu1(cosmo, filename, kgrid, TFnu_mass[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_nu2(cosmo, filename,  kgrid, TFnu2[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_extra(cosmo, filename, kgrid, TFextra[jz]);
        do_check(checknu_mass);
        checkaxion=gettransfer_axion(cosmo, filename, kgrid, TFaxion[jz]);
        do_check(checkaxion);
        fclose(fp);
    
    }

    free(filename);

    return 1;
}



int save_transfers_class(
    Cosmology *cosmo, 
    double zlist_transfer[], 
    double *kgrid,
    double **TFm, 
    double **TFgamma, 
    double **TFnu_massless, 
    double **TFnu1, 
    double **TFnu2, 
    double **TFextra
    ){
    //we will save the TF from CLASS

    int lengthname=200, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;
    int checkm, checknu, checknu_mass, checkgamma;

    for(jz=0;jz<Nz_transfer;jz++){
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/z%d_tk.dat", 
            boltzmann_tag, 
            cosmo->file_tag, 
            Nz_transfer-jz
        );
        //names start at z1=z_i, up to Nz_transfer-1 = z0. 
        //zlist transfer has been reversed.
        fp=fopen(filename, "r");
        checkm=gettransfer_matter(cosmo, filename, kgrid, TFm[jz]);
        do_check(checkm);
        checkgamma=gettransfer_gamma(cosmo, filename, kgrid, TFgamma[jz]);
        do_check(checkgamma);
        checknu=gettransfer_nu_massless(cosmo, filename, kgrid, TFnu_massless[jz]);
        do_check(checknu);
        checknu_mass=gettransfer_nu1(cosmo, filename, kgrid, TFnu1[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_nu2(cosmo, filename,  kgrid, TFnu2[jz]);
        do_check(checknu_mass);
        checknu_mass=gettransfer_extra(cosmo, filename, kgrid, TFextra[jz]);
        do_check(checknu_mass);
        fclose(fp);
    }

    free(filename);

    return 1;
}


////////////////////////////////////////////////////////////////////////////
//// These functions run Boltzmann over the redshift list zlist_transfer and 
//// saves transfer functions.  
////////////////////////////////////////////////////////////////////////////
int run_camb(Cosmology *cosmo, double zlist_transfer[]){

    //We will copy a version of params_base.ini and modify it
    system("cp CAMB_Current/params_base.ini params_collapse.ini");

    //create appropriate folder for CAMB output
    int lengthfilename=200;
    char *folder_name;
    folder_name=(char *)malloc((lengthfilename+1)*sizeof(char));

    //lengthfilename=sprintf(folder_name,"transfer_files_%d/test", file_tag);
    //fp=fopen(folder_name, "w");
    //fclose(fp);

    lengthfilename=sprintf(folder_name,"mkdir Boltzmann_%d", boltzmann_tag);
    system(folder_name);
    lengthfilename=sprintf(
        folder_name,
        "mkdir Boltzmann_%d/transfer_files_%d", 
        boltzmann_tag, 
        cosmo->file_tag
    );

    system(folder_name);

    int lengthname=200, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;

    lengthname=sprintf(filename,"params_collapse.ini");
    fp=fopen(filename, "a+");

    fprintf(fp, "ombh2 = %.4f \n",cosmo->omegab);
    fprintf(fp, "omch2 = %.4f \n",cosmo->omegac);
    //fprintf(fp, "omnuh2 = %.4f \n",cosmo->omeganu1);
    fprintf(fp, "omnuh2 = %.4f \n",3.*cosmo->omeganu1);
    fprintf(fp, "hubble = %.4f \n",100.*cosmo->h);
    //fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
    if(cosmo->omeganu1==0.){
        fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
        fprintf(fp, "nu_mass_eigenstates = %d \n", 1);
        fprintf(fp, "massive_neutrinos  = %d \n", 0);
    }
    else{
        fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
        fprintf(fp, "nu_mass_eigenstates = %d \n", 1);
        fprintf(fp, "massive_neutrinos  = %d \n", 3);
    }

    fprintf(fp, "scalar_amp(1) = %.2le \n", cosmo->As);
    fprintf(fp, "scalar_spectral_index(1) = %.4f \n", cosmo->ns);
    fprintf(fp, "transfer_num_redshifts = %ld \n", (long) Nz_transfer);
    for(jz=0;jz<Nz_transfer;jz++){
        fprintf(
            fp, 
            "transfer_filename(%d) = transfer_out_z%.3f \n", 
            jz+1, 
            zlist_transfer[jz]
        );
        fprintf(
            fp, 
            "transfer_redshift(%d) = %.3f \n", 
            jz+1, 
            zlist_transfer[jz]
        );
    }
    fprintf(
        fp, 
        "output_root = Boltzmann_%d/transfer_files_%d/ \n", 
        boltzmann_tag, 
        cosmo->file_tag
    );

    //fprintf(fp, "nu_mass_eigenstates = %d \n", cosmo->counter_massive_nus);
    //fprintf(fp, "massive_neutrinos  = %d \n", cosmo->counter_massive_nus);

    fclose(fp);

    //system("make clean");
    //system("make");

    system("mv params_collapse.ini CAMB_Current/params_collapse.ini");
    system("CAMB_Current/./camb CAMB_Current/params_collapse.ini\
         > output_cmd.txt");
    system("rm output_cmd.txt");

    free(filename);

    return 1;
}


int run_axioncamb(Cosmology *cosmo, double zlist_transfer[]){

    //We will copy a version of params_base.ini and modify it
    system("cp axionCAMB_Current/params_base.ini params_collapse.ini");

    //create appropriate folder for CAMB output
    int lengthfilename=200;
    char *folder_name;
    folder_name=(char *)malloc((lengthfilename+1)*sizeof(char));

    //lengthfilename=sprintf(folder_name,"transfer_files_%d/test", file_tag);
    //fp=fopen(folder_name, "w");
    //fclose(fp);

    lengthfilename=sprintf(folder_name,"mkdir Boltzmann_%d", boltzmann_tag);
    system(folder_name);
    lengthfilename=sprintf(
        folder_name,
        "mkdir Boltzmann_%d/transfer_files_%d",
        boltzmann_tag,
        cosmo->file_tag
    );

    system(folder_name);

    int lengthname=200, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;

    lengthname=sprintf(filename,"params_collapse.ini");
    fp=fopen(filename, "a+");

    fprintf(fp, "ombh2 = %.4f \n",cosmo->omegab);
    fprintf(fp, "omch2 = %.4f \n",cosmo->omegac);
    //fprintf(fp, "omnuh2 = %.4f \n",cosmo->omeganu1);
    fprintf(fp, "omnuh2 = %.4f \n",3.*cosmo->omeganu1);
    fprintf(fp, "hubble = %.4f \n",100.*cosmo->h);
    fprintf(fp, "omaxh2 = %.4le \n",cosmo->omega_ax);
    fprintf(fp, "m_ax = %.4le \n",cosmo->m_ax);
    //fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
    if(cosmo->omeganu1==0.){
        fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
        fprintf(fp, "nu_mass_eigenstates = %d \n", 1);
        fprintf(fp, "massive_neutrinos  = %d \n", 0);
    }
    else{
        fprintf(fp, "massless_neutrinos = %.4f \n", cosmo->Neff);
        fprintf(fp, "nu_mass_eigenstates = %d \n", 1);
        fprintf(fp, "massive_neutrinos  = %d \n", 3);
    }

    fprintf(fp, "scalar_amp(1) = %.2le \n", cosmo->As);
    fprintf(fp, "scalar_spectral_index(1) = %.4f \n", cosmo->ns);
    fprintf(fp, "transfer_num_redshifts = %ld \n", (long) Nz_transfer);
    for(jz=0;jz<Nz_transfer;jz++){
        fprintf(fp, "transfer_filename(%d) = transfer_out_z%.3f \n", jz+1, zlist_transfer[jz]);
        fprintf(fp, "transfer_redshift(%d) = %.3f \n", jz+1, zlist_transfer[jz]);
    }
    fprintf(fp, "output_root = Boltzmann_%d/transfer_files_%d/ \n", boltzmann_tag, cosmo->file_tag);

    //fprintf(fp, "nu_mass_eigenstates = %d \n", cosmo->counter_massive_nus);
    //fprintf(fp, "massive_neutrinos  = %d \n", cosmo->counter_massive_nus);

    fclose(fp);

    //system("make clean");
    //system("make");

    system("mv params_collapse.ini axionCAMB_Current/params_collapse.ini");
    system("axionCAMB_Current/./camb axionCAMB_Current/params_collapse.ini > output_cmd.txt");
    system("rm output_cmd.txt");

    free(filename);

    return 1;
}


int run_class(Cosmology *cosmo, double zlist_transfer[]){

    //We will copy a version of params_base.ini and modify it
    system("cp CLASS_Current/explanatory_base.ini explanatory_collapse.ini");

    //create appropriate folder for CAMB output
    int lengthfilename=300;
    char *folder_name;
    folder_name=(char *)malloc((lengthfilename+1)*sizeof(char));

    //lengthfilename=sprintf(folder_name,"transfer_files_%d/test", file_tag);
    //fp=fopen(folder_name, "w");
    //fclose(fp);

    lengthfilename=sprintf(folder_name,"mkdir Boltzmann_%d", boltzmann_tag);
    system(folder_name);
    lengthfilename=sprintf(
        folder_name,
        "mkdir Boltzmann_%d/transfer_files_%d", 
        boltzmann_tag, 
        cosmo->file_tag
    );

    system(folder_name);

    int lengthname=300, jz;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;

    lengthname=sprintf(filename,"explanatory_collapse.ini");
    fp=fopen(filename, "a+");

    fprintf(fp, "omega_b = %.4f \n",cosmo->omegab);
    fprintf(fp, "omega_cdm = %.4f \n",cosmo->omegac);
    fprintf(fp, "N_ur = %.4f \n",cosmo->Neff);
    fprintf(fp, "h = %.4f \n", cosmo->h);
    fprintf(fp, "T_cmb = %.4f \n", cosmo->T0_photons);
    fprintf(fp, "A_s = %.4le \n", cosmo->As);
    fprintf(fp, "n_s = %.4f \n", cosmo->ns);

    fprintf(fp, "z_pk = ");
    for(jz=0;jz<Nz_transfer-1;jz++){
        fprintf(fp, "%.3f, ", zlist_transfer[jz]);
    }

    fprintf(fp, "%.3f \n", zlist_transfer[Nz_transfer-1]);
    //last one apart, to not have a comma

    fprintf(
        fp, 
        "root = Boltzmann_%d/transfer_files_%d/ \n", 
        boltzmann_tag, 
        cosmo->file_tag
    );

    fprintf(
        fp, 
        "N_ncdm = %d \n", 
        cosmo->counter_massive_nus + (cosmo->Omega_extra>0)
    );

    if(cosmo->mnu1>0 || cosmo->omega_extra>0) {
        fprintf(fp, "m_ncdm = ");
    }
    if(cosmo->mnu1>0) {
        fprintf(fp, " %.3le ", cosmo->mnu1);
    }
    if(cosmo->mnu2>0) {
        fprintf(fp, " ,%.3le ", cosmo->mnu2);
    }
    if(cosmo->Omega_extra>0) {
        if(cosmo->mnu1>0) {
            fprintf(fp, " , ");
        }
        fprintf(fp, " %.3le ", cosmo->m_extra);
    }
    fprintf(fp, " \n");
    if(cosmo->mnu1>0 || cosmo->Omega_extra>0) {
        fprintf(fp, "T_ncdm = ");
    }
    if(cosmo->mnu1>0) {
        fprintf(fp, " %.3le ", cosmo->T0_nu/cosmo->T0_photons);
    }
    if(cosmo->mnu2>0) {
        fprintf(fp, " ,%.3le ", cosmo->T0_nu/cosmo->T0_photons);
    }
    if(cosmo->Omega_extra>0) {
        if(cosmo->mnu1>0) {
            fprintf(fp, " , ");
        }
        fprintf(fp, " %.3le ", cosmo->T0_extra/cosmo->T0_photons);
    }
    fprintf(fp, " \n");

    fclose(fp);

    //system("make clean");
    //system("make");

    system("mv explanatory_collapse.ini \
        CLASS_Current/explanatory_collapse.ini");
    lengthfilename=sprintf(
        folder_name,
        "CLASS_Current/./class CLASS_Current/explanatory_collapse.ini > Boltzmann_%d/transfer_files_%d/output_cmd.txt", 
        boltzmann_tag, 
        cosmo->file_tag
    );

    system(folder_name);

    //system("CLASS_Current/./class CLASS_Current/explanatory_collapse.ini > output_cmd.txt");
    //system("rm output_cmd.txt");

    free(filename);

    return 1;
}


////////////////////////////////////////////////////////////////////////////////
////    we RUN Boltzmann code to get the transfer functions between z=zi 
////    and z=0.    
////////////////////////////////////////////////////////////////////////////////

int boltzmann(Cosmology *cosmo, double *zlist_transfer){
    //this function runs Boltzmann code (CLASS, CAMB, AXIONCAMB), given z_collapse_array as
    //input, and stores in j_collapse_array the numbers of those z_collapses 
    //(since CLASS does not save z in output) and saves the list of zs in 
    //zlist_transfer.

    double *zlist_transfer_code;
    //same as zlist_transfer but reversed so it goes upwards. Because of CAMB...
    zlist_transfer_code=allocate_1D_array(Nz_transfer);

    int lengthname=200;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;

    long jz;

    const double zstep_sqrt_transfer=(sqrt(zf)-sqrt(zi))/(Nz_transfer-1); 
    //for reading the transfer functions
    //we use sqrt(z) spacing, since it samples better low z, and we cannot do 
    //log-spacing since we want to get to z=0.

    for(jz=0;jz<Nz_transfer;jz++){
        zlist_transfer_code[jz] = pow(sqrt(zi)+zstep_sqrt_transfer*jz,2.);
        //printf("%le \n",zlist_transfer_code[jz]);
    }

    //we tweak the zlist a bit
    //first we make sure the first z after zi is close to have accurate 
    //derivative
    double epsilon_z=0.005; 
    // dz/z at z_i to compute derivative, 0.005 is good, not terribly important,
    //since it's a small correction.

    zlist_transfer_code[0]=zi*(1.+epsilon_z); 
    //we overwrite the first non-zi z, for derivative at zi to be accurate.

    zlist_transfer_code[1]=zi;
    zlist_transfer_code[2]=zi*(1.-epsilon_z); 
    //we also overwrite the third non-zi z, for derivative at zi to be more 
    //accurate.

    //now we make sure that we have z_collapse exactly to avoid 
    //interpolation in z.
    cosmo->j_collapse_array[0]=find_value_reverse(
        Nz_transfer, 
        zlist_transfer_code, 
        cosmo->z_collapse_array[0]
    );

    //we do not want to remove z=0, so we check we are not:
    if(cosmo->j_collapse_array[0] == Nz_transfer-1){
        cosmo->j_collapse_array[0]--;
    }

    //we change whatever element of z was there for z_collapse.
    zlist_transfer_code[cosmo->j_collapse_array[0]]=cosmo->z_collapse_array[0];

    //and now we check that the z are not "colliding", and have 
    //different values.
    for(jz=1;jz<cosmo->N_zcoll;jz++){
        cosmo->j_collapse_array[jz]=find_value_reverse(
            Nz_transfer, 
            zlist_transfer_code, 
            cosmo->z_collapse_array[jz]
        ); 
        //printf("j=%ld, zlist=%.3le, and z_coll=%.3le \n",j_collapse, 
        //zlist_transfer_code[j_collapse],z_collapse);

        if(cosmo->j_collapse_array[jz]==cosmo->j_collapse_array[jz-1]){
            cosmo->j_collapse_array[jz]--;
        }
        //printf("z_coll=%f, z_code=%f, j_coll=%d \n\n\n",z_collapse_array[jz],
        //zlist_transfer_code[j_collapse_array[jz]],j_collapse_array[jz]);

        zlist_transfer_code[cosmo->j_collapse_array[jz]]=cosmo->z_collapse_array[jz];
    }

    for(jz=0;jz<Nz_transfer;jz++){
        zlist_transfer_code[jz] = (int)(zlist_transfer_code[jz] * 1000) / 1000.0f; 
        //we do this so the calculated values are exactly the ones with save, 
        //with the .3f
    }

    int boltzmann_check;

    if(boltzmann_tag == _CLASS_){//CLASS
        if(run_boltzmann_option!=0){//whether to run it
            boltzmann_check=run_class(cosmo, zlist_transfer_code);
        }
    }
    else if (boltzmann_tag == _CAMB_){//CAMB
        if(run_boltzmann_option!=0){//whether to run it
            boltzmann_check=run_camb(cosmo, zlist_transfer_code);
        }
    }
    else if (boltzmann_tag == _AXIONCAMB_){//AXIONCAMB
        if(run_boltzmann_option!=0){//whether to run it
            boltzmann_check=run_axioncamb(cosmo, zlist_transfer_code);
        }
    }
    else {
        printf("Choose either CLASS, CAMB, or axionCAMB \n");
        return -1;
    }

    do_check(boltzmann_check); //if 0 it has NOT run boltzmann

    //we reverse zlist_transfer now so it's ordered up 
    //(CAMB requires it to be ordered down...)

    for(jz=0;jz<Nz_transfer;jz++){
        zlist_transfer[jz]=zlist_transfer_code[Nz_transfer-1-jz];
        //printf("%le \n",zlist_transfer[jz]);
    }

    //we save the zlist to a file if we use CLASS.
    if(boltzmann_tag == _CLASS_){//CLASS
        lengthname=sprintf(
            filename,
            "Boltzmann_%d/transfer_files_%d/zlist.dat",
            boltzmann_tag,
            cosmo->file_tag
        ); //zlist for every tk

        fp=fopen(filename,"w");
        if(print_headers!=0){
            fprintf(fp, "z\t CLASS index j-1 \n");
        }
        for(jz=0;jz<Nz_transfer;jz++){
            fprintf(fp, "%le %ld \n",zlist_transfer_code[jz],jz);
        }
        fclose(fp);
    }

    free(zlist_transfer_code);

    return 1;
}


int read_input_file (char *filename, double parameter_values[Ninput]){
//this function reads the parameters from an input file in filename
//spits out parameter values.

    int i;

    std::string buffer;

    std::ifstream is(filename);    // open file

    char c;
    while (is.get(c)){    // loop getting single characters
        buffer += c;
    }

    is.close();
    // close file

    //to check it has read it:
    //std::cout << buffer << "\n";

    std::string parameter_lookup[] = {
            "N_Mhalo",
            "N_zcoll",
            "N_klong",
            "tag_thermal_relic",
            "tag_sterile_nu",
            "Mhalo_min",
            "Mhalo_max",
            "z_collapse_bot",
            "z_collapse_top",
            "kbot",
            "ktop",
            "m_TR",
            "omega_TR",
            "m_SN",
            "omegab",
            "omegac",
            "hubble",
            "A_s",
            "n_s",
            "mnu1",
            "mnu2",
            "N_eff_input",
            "omega_ax", 
            "m_ax" 
    };

    std::string parameter_names[Ninput] = {};
    for(i=0; i<Ninput; i++){
        parameter_names[i]=parameter_lookup[i];
    } 
    std::string param_name;
    double param_value;
    for (i=0;i<Ninput;i++){
        param_name = parameter_names[i];
        param_value = read_param(buffer, param_name);
        if(i<Ninput_int){    //if it's one of the int parameters
            parameter_values[i] = (int) param_value;
        }
        else{
            parameter_values[i] = param_value;
        }
        std::cout << param_name << " = " << param_value << "\n";
    }
    return 1;
}


double read_param (std::string buffer, std::string parameter_name){
//reads parameter value from buffer, after finding parameter_name

    double parameter_value;    //the output

    int i;

    int length_line = 50;
    char line[length_line];    //temporary variable, to store the line where the 
                               //variable name appears
    std::size_t position = buffer.find(parameter_name);    //we find
                                                           //in buffer

    for(i=0;i<length_line;i++){
        line[i]    = buffer.at(position+i);     //we save the parameter name and  
                                                //its value in temporary line
    }

    if (debug_mode > 1){
        std::cout << line << "\n";    //print it on screen just to see it
    }

    sscanf(line, "%*[^=]= %le ", &parameter_value);    //save the number after
                                                       //the equal sign
    //to check it's working:
    if (debug_mode > 0){
        std::cout << parameter_name;
        printf("= %le \n", parameter_value);
        std::cout << "\n";
    }
    return parameter_value;
}


int prepare_cosmology(Cosmology *cosmo, double *parameter_values){
//we assign the read parameters.

    cosmo->N_Mhalo = parameter_values[0];
    cosmo->N_zcoll = parameter_values[1];
    cosmo->N_klong = parameter_values[2];
    cosmo->tag_thermal_relic = parameter_values[3];
    cosmo->tag_sterile_nu = parameter_values[4];
    cosmo->Mhalo_min = parameter_values[5];
    cosmo->Mhalo_max = parameter_values[6];
    cosmo->z_collapse_bot = parameter_values[7];
    cosmo->z_collapse_top = parameter_values[8];
    cosmo->kbot = parameter_values[9];
    cosmo->ktop = parameter_values[10];
    cosmo->m_TR = parameter_values[11];
    cosmo->omega_TR = parameter_values[12];
    cosmo->m_SN = parameter_values[13];
    cosmo->omegab = parameter_values[14];
    cosmo->omegac = parameter_values[15];
    cosmo->h = parameter_values[16];
    cosmo->As = parameter_values[17];
    cosmo->ns = parameter_values[18];
    cosmo->mnu1 = parameter_values[19];
    cosmo->mnu2 = parameter_values[20];
    cosmo->Neff_input = parameter_values[21];

    if (boltzmann_tag==_AXIONCAMB_){ 
        cosmo->omega_ax = parameter_values[22]; 
        cosmo->m_ax = parameter_values[23]; 
    }

    cosmo->counter_massive_nus = (cosmo->mnu1>0) + (cosmo->mnu2>0); 
    //we subtract one from Neff per neutrino, and we use it for 
    //counting columns in CLASS

    cosmo->Neff = (
        cosmo->Neff_input - constant_change_Neff * (
            cosmo->counter_massive_nus 
            + cosmo->tag_sterile_nu*(cosmo->m_SN<m_SN_thresh_Neff) 
        )
    ); 
    //we also subtract 1 if SN is below thresh mass, as we assume it's a 3rd
    //regular neutrino

    cosmo->file_tag = (
        (cosmo->tag_thermal_relic>0) 
        * (int)(1000.*cosmo->omega_TR + 1000*cosmo->m_TR + 0.5f)
        + (cosmo->tag_sterile_nu>0) * (int)(100*cosmo->m_SN + 0.5f)
    );
    cosmo->file_tag = (0);  
    //some placeholder for what to use as filetag. + 0.5f so it rounds

    cosmo->omega_SN = cosmo->m_SN/mass_constant_nu; 
    //Omega h^2, since Temperature is the same than regular nus.

    cosmo->omega_extra = (
        (cosmo->tag_thermal_relic>0)*cosmo->omega_TR 
        + (cosmo->tag_sterile_nu>0)*cosmo->omega_SN
    );

    cosmo->Omega_extra = cosmo->omega_extra/cosmo->h/cosmo->h;

    cosmo->m_extra = (
        (cosmo->tag_thermal_relic>0)*cosmo->m_TR 
        + (cosmo->tag_sterile_nu>0)*cosmo->m_SN
    );

    cosmo->omeganu1=cosmo->mnu1/mass_constant_nu; //for mnu1 in eV.

    cosmo->omeganu2=cosmo->mnu2/mass_constant_nu; //for mnu2 in eV.

    cosmo->Omeganu1=cosmo->omeganu1/cosmo->h/cosmo->h;

    cosmo->Omeganu2=cosmo->omeganu2/cosmo->h/cosmo->h;

    cosmo->Omega_ax = cosmo->omega_ax/cosmo->h/cosmo->h; //for axion

    cosmo->OmegaM=(cosmo->omegab+cosmo->omegac)/cosmo->h/cosmo->h; 
    //b+c energy density, this is what haloes are made of (and what enters R_M)

    cosmo->T0_photons = 2.7255; 
    //temperature of photons today, in K, from Planck2015 (I assume this is not
    //going to be changed often)

    cosmo->T0_nu = cosmo->T0_photons * ratio_T_nu_T_gamma;
    //temperature of neutrinos today, in K.

    cosmo->T0_TR = (
        (cosmo->tag_thermal_relic>0) ? cosmo->T0_nu 
        * pow(cosmo->omega_TR * mass_constant_nu/cosmo->m_TR, 1./3) : 0.0
    ); 
    //temperature of TR today in K. //Fixed through cosmic abundance,
    //assuming relativistic at decoupling and non-rel today

    cosmo->T0_extra = (
        (cosmo->tag_sterile_nu>0)*cosmo->T0_nu 
        + (cosmo->tag_thermal_relic>0)*cosmo->T0_TR
    ); 
    //SN have nu temperature.

    cosmo->OmegaG = (
        omega_constant/cosmo->h/cosmo->h 
        * pow(cosmo->T0_photons,4.0)
    );
    //Omega_photons at z=0

    cosmo->Omeganu_massless = (
        cosmo->OmegaG * 7./8 * cosmo->Neff 
        * pow(cosmo->T0_nu/cosmo->T0_photons,4.0)
    ); 
    //Omega_masslessneutrinos

    cosmo->OmegaR = cosmo->OmegaG + cosmo->Omeganu_massless; 
    //radiation (at all z)

    cosmo->OmegaL = (
        1.0 - cosmo->OmegaM - cosmo->Omeganu1 - cosmo->Omeganu2 
        - cosmo->OmegaR - cosmo->Omega_extra - cosmo->Omega_ax
    ); 
    //close the Friedmann eq., no curvature.

    cosmo->H0_Mpc=100.*cosmo->h/(c_light); 
    //H0 in Mpc-1

    //we save the zs and masses we calculate for
    //masses:
    cosmo->Mhalo_array = allocate_1D_array(cosmo->N_Mhalo);
    int iM;
    double logstep_Mhalo=0;//if only one mass take the lowest in the interval.

    if(cosmo->N_Mhalo>1){
        logstep_Mhalo=log(cosmo->Mhalo_max/cosmo->Mhalo_min)/(cosmo->N_Mhalo-1.);
    }

    for(iM=0;iM<cosmo->N_Mhalo;iM++){//we populate the z_collapse_array
        cosmo->Mhalo_array[iM] = cosmo->Mhalo_min * exp(logstep_Mhalo*iM);
        if (debug_mode>0) printf("Mhalo[%d]=%.1le Msun \n",iM,cosmo->Mhalo_array[iM]);
    }

    //redshifts:
    cosmo->j_collapse_array = allocate_1D_array_int(cosmo->N_zcoll); 
    //these are the indices, we find them with boltzmann()
    
    cosmo->z_collapse_array = allocate_1D_array(cosmo->N_zcoll);
    int iz;
    double z_collapse_step=0.;
    if(cosmo->N_zcoll>1){
        z_collapse_step=(
            (cosmo->z_collapse_top-cosmo->z_collapse_bot)
            / (cosmo->N_zcoll-1.)
        );
    }
    for(iz=0;iz<cosmo->N_zcoll;iz++){//we populate the z_collapse_array
        cosmo->z_collapse_array[iz] = cosmo->z_collapse_bot + z_collapse_step*iz;
    }
    //note that we do not initialize z_collapse.

    //and also how many k inputs:
    cosmo->klong_list_input = allocate_1D_array(cosmo->N_klong);

    int ik;
    double logstep_klong = 0.0;

    if((cosmo->N_klong) > 1){
        logstep_klong=log(cosmo->ktop/cosmo->kbot)/(cosmo->N_klong-1.);
    }
    for(ik=0;ik<cosmo->N_klong;ik++){ //we populate the z_collapse_array
        cosmo->klong_list_input[ik] = cosmo->kbot * exp(logstep_klong*ik);
    }
    return 1;
}
